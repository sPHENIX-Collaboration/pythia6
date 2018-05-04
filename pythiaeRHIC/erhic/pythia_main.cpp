#include <iostream>
#include <sstream>
#include <string>

// GNU long option variation of getopt, for command line option parsing.
#include <getopt.h>

// ROOT headers.
#include <TFile.h>
#include <TObjString.h>
#include <TTree.h>

#include "factory.h"
#include "pythia_commons.h"
#include "pythia_erhic.h"

// Forward declarations
std::string initialise(int argc, char* argv[]);
void write_run_data(TFile&);

/**
 C++ wrapper for calling the eRHIC PYTHIA 6 implementation
 and building eRHIX C++ events to store in a ROOT tree.
 */
int main(int argc, char* argv[]) {
   std::string filename = initialise(argc, argv);

   // Open the output file/tree and create an event factory,
   // which we use to add the event branch to the tree.
   TFile file(filename.c_str(), "recreate");
   TTree tree("EICTree", "Output direct from pythia!");
   Factory factory;
   factory.Branch(tree, "event");

   std::stringstream stream;
   // The PYTHIA6 generate() subroutine handles tracking the event
   // count. Continue looping until it returns a non-zero value,
   // indicating the requested number of events have been generated.
   while(__pythia6_MOD_generate() == 0) {
      erhic::EventPythia* event = factory.Create();
      event->SetN(tree.GetEntries() + 1);
      tree.Fill();
      // \todo Add radiative photon
      std::cout << std::endl;
   } // for
   file.Write(); // Writes the ROOT TTree
   write_run_data(file);
   __pythia6_MOD_finish();
   return 0;
}

/**
 Initialise the program, processing command line arguments.
 Returns the name of the ROOT file to generate.
 */
std::string initialise(int argc, char* argv[]) {
   std::string filename("pythia.root");
   // Parse command line options.
   __pythia6_MOD_printascii = 1; // Generate ASCII output by default
   static struct option long_options[] = {
      // These set a flag.
      {"ascii", no_argument, &__pythia6_MOD_printascii, 1},
      {"noascii", no_argument, &__pythia6_MOD_printascii, 0},
      // These don't set a flag.
      // The arguments are processed below.
      {"out", required_argument, NULL, 'o'},
      {NULL, 0, NULL, 0}
   };
   // Loop through options.
   int option_index = 0;
   int code(0);
   while((code = getopt_long(argc, argv, "o:r:",
                             long_options, &option_index)) not_eq -1) {
      switch(code) {
         case 0:
            if(long_options[option_index].flag not_eq 0) {
               break;
            } // if
            printf("option %s", long_options[option_index].name);
            if(optarg) {
               printf (" with arg %s", optarg);
            } // if
            printf("\n");
            break;
         case 'o':
            filename = optarg;
            break;
         default:
            abort();
      } // switch
   } // while
   // Now that we have processed all command line arguments, call
   // the Fortran PYTHIA6 initialisation routine.
   // This handles processing options from the steer file.
   // If the return value is non-zero, exit.
   int result = __pythia6_MOD_initialise();
   if(0 not_eq result) {
      exit(result);
   } // if
   return filename;
}

/**
 Extract the following run data from PYTHIA and write as
 TObjStrings to the provided file with the names shown:
 <ul>
   <li>Total sampled cross section as "crossSection"</li>
   <li>Total number of generated events as "nEvents"</li>
   <li>Total number of event trials as "nTrials"</li>
 </ul>
*/
void write_run_data(TFile& file) {
   TObjString text;
   std::stringstream stream;
   // Total cross section is stored in PYTHIA in pari(1) in millibarns.
   // Multiply by 1000 to get it in microbarn.
   stream << pari(1) * 1000.;
   text.SetString(stream.str().c_str());
   file.WriteObject(&text, "crossSection");
   // Total number of generated events is stored in PYTHIA in msti(5).
   stream.str("");
   stream.clear();
   stream << msti(5);
   text.SetString(stream.str().c_str());
   file.WriteObject(&text, "nEvents");
   // Total number of generateds is stored in PYTHIA in ngen(0, 3).
   stream.str("");
   stream.clear();
   stream << ngen(0, 3);
   text.SetString(stream.str().c_str());
   file.WriteObject(&text, "nTrials");
}

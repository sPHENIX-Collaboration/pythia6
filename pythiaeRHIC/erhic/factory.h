/**
 \file factory.h
 
 \author Thomas Burton 
 \date 10/9/12
 \copyright 2012 BNL.
 */

#ifndef _PYTHIAERHIC_FACTORY_H_
#define _PYTHIAERHIC_FACTORY_H_

#include <string>

#include <eicsmear/erhic/EventFactory.h>

class TBranch;
class TTree;

// Forward declaration of erhic::EventPythia.
namespace erhic {
   class EventPythia;
} // namespace erhic

/**
 Factor class for PYTHIA events.
 Handles creating the ROOT tree branch and copying values from
 PYTHIA common blocks/modules to the event.
 */
class Factory : public erhic::VirtualEventFactory {
public:
   
   /** Constructor */
   Factory();
   
   /** Destructor */
   virtual ~Factory();
   
   /** Populate the stored event with the current contents of
    PYTHIA and return a pointer to the event.
    Do not delete this pointer!
   */
   erhic::EventPythia* Create();
   
   /**
    Returns a string with the full (including namespace) class name
    of the event type produced.
    This is important for use with ROOT TTree to ensure the correct
    event type in branches.
    */
   virtual std::string EventName() const;
   
   /**
    Add a branch named "name" for the event type generated
    by this factory to a ROOT TTree.
    Returns a pointer to the branch, or NULL in the case of an error.
    */
   virtual TBranch* Branch(TTree& tree, const std::string& name);
   
private:

   erhic::EventPythia* mEvent;
}; // class Factory

#endif // _PYTHIAERHIC_FACTORY_H_

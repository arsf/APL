//------------------------------------------------------------------------- 
//Copyright (c) 2013 Natural Environment Research Council (NERC) UK 
// 
//This file is part of APL (Airborne Processing Library)
//Licensed under the APL Open Software License version 1.0 
// 
//You should have received a copy of the Licence along with the APL source 
//If not, please contact arsf-processing@pml.ac.uk 
//-------------------------------------------------------------------------

/**************************************************************************
*     Class to handle command lines, sort options/arguments. 
***************************************************************************/

#include "commandline.h"

#ifndef COMMANDLINEDEBUG
   #define DEBUGPRINT(X) 
#else
   #define DEBUGPRINT(X) std::cout<< X <<std::endl;
#endif

//Conctructor for commandline object. Pass it the argv[] and argc parameters
//fills the elements map with the options and arguements (by calling Handle function)
CommandLine::CommandLine(char** args,int numargs)
{
   //DEBUG output
   DEBUGPRINT("Entering CommandLine Constructor...");
   clgood=true;//set this to true initially
   argcount=NULL; //assign pointers to NULL
   bool retval=Handle(args,numargs);//handle the arguments and options
   if(retval==0)
   {
      //A problem occurred
      clgood=false;//set flag
      throw CommandLineException();//throw an exception
   }
   
   CountArgs();
}

//Destructor
CommandLine::~CommandLine()
{
   //DEBUG output
   DEBUGPRINT("Entering CommandLine destructor...");
   if(argcount!=NULL)
      delete[] argcount;
   
   //this->elements.clear();
}

//function to check command line options against avopts (strings of available options eg "-f")
//returns -ve number if options on cl are not in avopts [value == to number of incorrect options ] 
//returns 0 otherwise
int CommandLine::CheckAvailableOptions(const std::string* avopts,const int numstrings,std::string* const badopts)
{
   //DEBUG output
   DEBUGPRINT("Checking command line's available options...");
   std::map<std::string,std::string>::iterator it;
   std::string opt("");
   int count=0;
   bool available=false;
   //for each command line option check it exists in strings avopts
   for(it=elements.begin();it!=elements.end();it++)//for each cl option
   {
      opt=it->first;//get option from command line 
      for(int i=0;i<numstrings;i++)//for each possible cl option
      {
         if(!opt.compare(avopts[i]))
         {
            //OPtion is in available options
            //dont need to do anything but break out of loop
            available=true;
            break;
         }
         else
         {
            //Option is not in available options
            available=false;
         }
      }      
      if(available==false)
      {
         //increase count of bad options
         if(badopts!=NULL)
         {
            if((*badopts).find(opt)==std::string::npos)
               *badopts+=opt+" ";
         }
         count++;
      }
   }
   return -count;
}


bool CommandLine::IsOpt(const std::string str)
{
   int index=str.find_first_of('-');
   if(index==0)
   {
      //could be an option or -ve number
      std::string nums="-+0123456789.";
      if(str.find_first_not_of(nums)==std::string::npos)
         return false;  //only chars from nums - assume it is a number
      else
         return true; //assume an option
   }
   else
   {
      //not an option
      return false;
   }
}

//Function to iterate through the numargs of argv[] and sort into options and arguements
//Adding functionality so that -ve numbers can be used as args but not as opts
bool CommandLine::Handle(char** args,int numargs)
{
   //DEBUG output
   DEBUGPRINT("Entering CommandLine Handler function...");
   //for each element in args we need to ascertain whether it is a option or argument
   //define options as starting with '-' anything else is an argument
   std::string temparg;
   std::string tempargp1;
   std::string key(""),value("");
   //int index=0;

   //Add the executable name to the list
   this->strExecutableName=args[0]; 

   for(int arg=1;arg<numargs;arg++)
   {
      temparg=args[arg];
      //Add the argument to the formatted command line string 
      this->strFormattedCL+=" "+temparg;

      //index=temparg.find_first_of('-');
      //if(index==0)
      if(this->IsOpt(temparg))
      {
         //This is an option
         //Check if appending values
         if(key.compare(""))
         {
            elements[key]=value;
            //reset key and value
            key="";
            value="";
         }

         //Need to check next arg to see if that is an option or argument
         if(arg<numargs-1)
         {
            //Then its not the last argument so safe to check next one
            tempargp1=args[arg+1];
            //if(tempargp1.find_first_of('-')==0)
            if(this->IsOpt(tempargp1))
            {
               //Next is an option therefor no matching argument
               elements[temparg]=optiononly;
            }
            else
            {
               //next is an argument
               key=temparg;//set the key as the option
            }
         }
         else
         {
            //This is the last argument of command line
            //Enter this as option only (no matching arguments)
            elements[temparg]=optiononly;
         }
      }
      else
      {
         //This is not an option  
         if(!key.compare(""))
         {
            //this is a "lone" argument with no option ... this is not allowed
            std::cout<<"Lone argument. All arguments should match with an option."<<std::endl;
            return 0;
         }
         else
         {
            //Add string to the value
            if(!value.compare(""))//value is empty
               value=temparg;
            else //value is not empty
               value=value+';'+temparg;
   
            if(arg==(numargs-1)) //last one so need to add here (because this is a baaad loop . should fix it
               elements[key]=value;
         }
      }
   }
   return 1;
}

//function to count the number of arguments per option
void CommandLine::CountArgs()
{
   //This should only be run once
   if(argcount!=NULL)
   {
      //Has already been run..dont want to re-run it
      return;
   }
   //DEBUG output
   DEBUGPRINT("Counting commandline arguments...");
   std::string argstr="";
   std::map<std::string,std::string>::iterator it;
   int totalopts=this->elements.size(); // number of options
   if(totalopts==0) 
      return; //no options

   argcount=new int[totalopts]; //assign array of new ints
   int opt=0;//counter for the number of options
   for(it=this->elements.begin();it!=this->elements.end();it++)
   {
      //Get the arguments for the option
      argstr=it->second;
      if(!argstr.compare(optiononly))   
         argcount[opt]=0; //no argument just a dummy value
      else
      {
         //Count the numebr of args separated by ';'
         argcount[opt]=TotalOccurence(argstr,';') + 1; //+1 since there is one less ';' than argument
      }
  
      opt++;//increase opt
   }

}

// returns the number of arguments for the given key
int CommandLine::NumArgsOfOpt(const std::string key)
{
   std::map<std::string,std::string>::iterator it;
   int opt=0;//counter for the number of options
   std::string elementkey; // key from map
   if(argcount==NULL)
   {
      //throw CommandLineException("Cannot access number of arguments...argcount is NULL pointer");
      return -1;
   }
   for(it=this->elements.begin();it!=this->elements.end();it++)
   {
      elementkey=it->first;
      //Chec if matches key
      if(!key.compare(elementkey))
      {
         //key matches element from map   
         return argcount[opt];
      }
      opt++;//increase opt
   }  
   //No match for key in map
   return -1; 
}

//Function to return the argument when more than one for the option
std::string CommandLine::GetArg(const std::string opt,int argn)
{
   //DEBUG output
   DEBUGPRINT("Retrieving argument "<<argn<<" from value relating to option "<<opt);
   std::map<std::string,std::string>::iterator it;
   it=this->elements.find(opt);
   if(it==this->elements.end()) //no key opt found in map return emptty string
   {
      //throw CommandLineException("Error no option "+opt+" on command line");
      return "";
   }
   else //key found
   {
      if(argn>=NumArgsOfOpt(opt))      
      {
         return ""; //trying to access an argument numbered higher than the number of arguments
      }
      else
      {
         //Need to loop to get the argn-th argument         
         std::string strvalue=it->second;//get the value relating to the key
         unsigned int index=strvalue.find_first_of(';'),previndex=0;
         //if we want the first argument then can return it now
         if(argn==0)
            return strvalue.substr(previndex,index-previndex);//return the argument substring
         //else we want to keep looping until we get the required argument
         for(int i=0;i<argn;i++)//for each arg after 0 upto the one we want
         {
            previndex=index;
            index=strvalue.find(';',index+1);//get the index of the ';'
            if(index==std::string::npos)
            {
               if(i<(argn-1))
                  throw CommandLineException("Trying to retrieve a command argument greater than number of items present. Wanted item "+ToString(argn+1)+" of "+ToString(i));
               else
                  return strvalue.substr(previndex+1);
                  //index=strvalue.length()-1;
            }
         }
         return strvalue.substr(previndex+1,index-(previndex+1));//return the substring that = argument
      }
   }
}

//Function to review the command line options. Takes a stringstream to report to and
//array of available options, and size of the array
std::string CommandLine::ReviewCL(const std::string* availopts,const int number_of_possible_options,bool printout)
{
   std::stringstream info;
   //Ouput the number of options on command line
   if(printout==true)
   {
      info<<"There are "<<this->elements.size()<<" options on the command line."<<std::endl;
      info<<"Listing all available options:"<<std::endl;
   }

   int numargs=0;
   std::map<std::string,std::string>::iterator it;
   //For each available option list if it is present and its arguments
   for(int i=0;i<number_of_possible_options;i++)
   {
      numargs=this->NumArgsOfOpt(availopts[i]);
      if(numargs==-1)
      {
         if(printout==true)
         {
            info<<availopts[i]<<"\t\t"<<"This option is not on command line."<<std::endl;  
         }
      }
      else
      {
         if(printout==true)
         {
            info<<availopts[i]<<"\t\t"<<numargs<<" arguments for this option present."<<std::endl;
         }
         //Output each argument for this option
         for(int arg=0;arg<this->NumArgsOfOpt(availopts[i]);arg++)
            info<<this->GetArg(availopts[i],arg)<<" ";
         info<<std::endl;
      }

   }

   int numbadopts=this->CheckAvailableOptions(availopts,number_of_possible_options);
   if(numbadopts!=0)
   {
      info<<"There are "<<(-numbadopts)<<" bad (unknown) options on the command line."<<std::endl; 
      //throw CommandLine::CommandLineException(info.str());
   }

   return info.str();
}

bool CommandLine::OnCommandLine(const std::string opt)//check if string opt is present on command line
{
   std::map<std::string,std::string>::iterator it;      
   it=this->elements.find(opt);
   if(it!=this->elements.end()) 
      return true;
   else
      return false;
}

std::string CommandLine::GetArg(const std::string opt)
{
   std::map<std::string,std::string>::iterator it;
   it=this->elements.find(opt);
   if(it!=this->elements.end())
      return it->second;
   else
   {
      //throw CommandLineException("Error no option "+opt+" on command line");
      return "";
   }
}  //Returns the value relating to the key (opt) in the elements map


//----------------------------------------------------------------------------
//Describe the usage of the program - i.e. a help document
//----------------------------------------------------------------------------
std::string CommandLine::ProgramUsage(const int number_of_possible_options,const std::string* availableopts,const std::string* optsdescription)
{
   std::stringstream buffer;
   #ifdef _W32
      const int maxsetlength=77;
   #else
      const int maxsetlength=100;
   #endif

   size_t chartowrite=0;
   //Do some pre-formatting 
   //Find the longest string in the arguments list
   int strLength=0,maxlength=0;
   for(int o=0;o<number_of_possible_options;o++)
   {
      strLength=availableopts[o].length();
      if(strLength > maxlength)
         maxlength=strLength;
   }   

   //Maximum length of a line after the "   optionname   | " bit
   int linemax=maxsetlength-maxlength-6;

   if(linemax < 0)
      linemax=maxsetlength; //just use the maximum length and suffer the poor formatting

   //Start by outputing a preamble about what the program does and when it was compiled
   buffer<<"   Usage for: "<<strExecutableName<<std::endl;
   //Now add the arguments and their description
   buffer<<std::endl<<"Arguments:"<<std::endl;
   for(int o=0;o<number_of_possible_options;o++)
   {
      buffer<<"   "<<availableopts[o];
      for(unsigned int i=0;i<maxlength-availableopts[o].length();i++)
         buffer<<" ";
      buffer<<" | ";
      chartowrite=optsdescription[o].length();
      size_t startchar=0,spacebreak=0;
      while(chartowrite>(unsigned)linemax)
      {
         spacebreak=optsdescription[o].substr(startchar,linemax-1).rfind(' ');
         for(size_t c=0;c<spacebreak;c++)
         {
            buffer<<optsdescription[o].at(startchar+c);
         }
         buffer<<std::endl<<"   ";
         for(int i=0;i<maxlength;i++)
            buffer<<" ";
         buffer<<" | ";
         chartowrite=chartowrite-(spacebreak+1);
         startchar=startchar+(spacebreak+1);
      }

      for(size_t c=0;c<chartowrite;c++)
      {
         buffer<<optsdescription[o].at(startchar+c);
      }
      buffer<<std::endl;
   
   }

   //return the buffer to the terminal/file
   return buffer.str();
}

//----------------------------------------------------------------------------
// Return the formatted command line string
//----------------------------------------------------------------------------
std::string CommandLine::ReturnCLAsString()
{
   return strExecutableName+strFormattedCL;
}

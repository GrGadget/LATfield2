/*! \file LATfield2_SettingsFile.hpp
 \brief LATfield2_SettingsFile.hpp contain the class SettingsFile definition.
 \author N. Bevis
 */ 

#include "LATfield2_SettingsFile.hpp"
#include "LATfield2_parallel2d.hpp"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

namespace LATfield2
{
using std::endl;
using std::cerr;
using std::cout;


//CONSTANTS===========================
int SettingsFile::noCreate = 1;
int SettingsFile::autoCreate = 0;


//CONSTRUCTORS========================
SettingsFile::SettingsFile() 
{

#ifndef SERIAL
  isRoot_=parallel.isRoot();
#else
  isRoot_=true;
#endif

}

SettingsFile::SettingsFile(const std::string filename, const int mode, const int argc, char** argv) 
{

#ifndef SERIAL
  isRoot_=parallel.isRoot();
#else
  isRoot_=true;
#endif

  this->open(filename, mode, argc, argv);
}

//DESTRUCTOR==========================
SettingsFile::~SettingsFile() {this->close();}

//OPEN================================
void SettingsFile::open(const std::string filename, const int mode, const int argc, char** argv)
{
  char c;

  filename_=filename;
  mode_=mode;

  if(isRoot_)
    {      
      //Open file  
      file_.open(filename_.c_str(), std::fstream::in);
      if(!file_.is_open())
	  {
		  if((mode_ & SettingsFile::noCreate) == 0)
		  {
			  std::cout<<"SettingsFile: "<<filename_<<" not found."<<std::endl;
			  std::cout<<"SettingsFile: Creating..."<<std::endl;
			  this->create(filename_);
			  std::cout<<"creating ok"<<std::endl;
		  }
		  else
		  {
			  std::cout<<"SettingsFile: "<<filename_<<" not found and auto create off."<<std::endl;
			  std::cout<<"SettingsFile: Exiting..." << std::endl;
#ifndef SERIAL
			  parallel.abortRequest();
#else
			  exit(555);
#endif
		  }
	  }	
  
      //Read command line into stringstream
      for(int i=0; i<argc; i++)
		{
			for(int j=0; argv[i][j]!='\0'; j++)
			{
				stream_<<argv[i][j];
			}
			stream_<<endl;    
		} 
  
      //Read file into stringstream
      while(!file_.eof())
	  {
		  c=file_.get();
		  if(c=='#')
		  {
			  while(!file_.eof() && c!='\n') { c=file_.get(); }    
			  if(file_.eof()) { break; }
		  }
		  stream_.put(c);       
	  }
	  file_.close();  
    } 


#ifndef SERIAL
  //Broadcast results to all processes 
	
	
	parallel.barrier();
	
  if(parallel.size()>1)
    {
      if(parallel.isRoot())   
	  {
		  int len = stream_.str().length();
		  char* streamString = new char[len+1];
		  for(int i=0;i<=len;i++) { streamString[i]=stream_.str()[i]; }
		  parallel.broadcast(len, parallel.root());
		  parallel.broadcast(streamString, len+1, parallel.root());
	  }    
      else
	  {
		  int len;
		  char* streamString;
		  parallel.broadcast(len, parallel.root());
		  streamString=new char[len+1];
		  parallel.broadcast(streamString, len+1, parallel.root());
		  stream_<<streamString;
	  }    
	}
    
#endif
}

//FILE CLOSE============================
void SettingsFile::close()
{
  if(isRoot_) { filename_="."; }
}

//FILE CREATE===========================
void SettingsFile::create(const std::string filename)
{
    
  if(isRoot_)
    {
      filename_=filename;
      mode_=autoCreate;
      
      file_.open(filename_.c_str(), std::fstream::out);
      if(!file_.is_open())
	  {
		  std::cout<<"SettingsFile: Cannot create: "<<filename<<std::endl;
		  std::cout<<"SettingsFile: Exiting..."<<std::endl;
#ifndef SERIAL
		  parallel.abortRequest();
#else
		  exit(555);
#endif	
	  }
      else
	  {
		  file_.close();
		  file_.clear();
		  file_.open(filename.c_str(), std::fstream::in);
	  }
    }

  //parallel.barrier();
}

//PARAMETER READ===========================
template<class TemplateClass>
void SettingsFile::read(const std::string parameterName, TemplateClass &parameter)
{
  if(this->search(parameterName+'='))
    {
      stream_>>parameter;
    }
  else
    { 
      if(isRoot_)
	  {
		  //verifiy that the parameter name is no autocreate
		  if(parameterName!="autocreate")
		  {
			  std::cout<<"SettingsFile: "<<filename_<<" has no parameter: "<<parameterName<<std::endl;
			  std::cout<<"SettingsFile: No command-line override given either"<<std::endl;  
			  
			  if((mode_ & SettingsFile::noCreate) == 0 )
			  {
				  std::cout << "SettingsFile: Adding with current value: " << parameter << std::endl; 
				  this->add(parameterName, parameter);
			  }
			  else
			  {
				  std::cout << "SettingsFile: Auto create off. Exiting..." << std::endl;
#ifndef SERIAL
				  parallel.abortRequest();
#else
				  exit(555);
#endif
			  }
		  }		  
		  
	  }
#ifndef SERIAL
		parallel.barrier();
#endif
    }
}

//PARAMETER WRITE===========================
template<class TemplateClass>
void SettingsFile::add(const std::string parameter_name, const TemplateClass &parameter)
{
  if(isRoot_)
    {
      file_.clear();
      file_.open(filename_.c_str(), std::ios::out | std::ios::app);
      file_ << parameter_name << '=' << parameter << std::endl;
      if(!file_.good())
	  {
	  std::cout << "SettingsFile: Could not write to file: " << filename_ << std::endl;
	  std::cout << "SettingsFile: Exiting... " << std::endl;
#ifndef SERIAL
		parallel.abortRequest();
#else
		exit(555);
#endif
	  }
      file_.close();
    }

}
template<class TemplateClass>
void SettingsFile::write(const std::string parameter_name, const TemplateClass &parameter)
{
	if(isRoot_)
    {
		unsigned int i=0;
		unsigned int l=0;
		int line_number=0;
		
		char c;
		std::string line_temp;
		
        if(file_.good())file_.close();
		file_.clear();
		file_.open(filename_.c_str(), std::ios::in);
		
		
		//get number of non void line
		file_.seekg(0);
		while (!file_.eof())
		{
			getline(file_,line_temp);
			line_number++;
		}
		
		//creat line array and this_param array
        std::string * lines;
        lines = new std::string[line_number];
        
		bool this_param[line_number];
		file_.seekg(0);
		file_.close();
		file_.open(filename_.c_str(), std::ios::in);
		
		//implement non void lines
		l=0;
		while (!file_.eof())
		{
			getline(file_,line_temp);
			
            lines[l]=line_temp;
            l++;
			
		}
		for(l=0;l<line_number;l++)
		{
			for(i=0;i<parameter_name.length();i++)
			{
				c = lines[l][i];
				if(c==parameter_name[i])this_param[l]=true;
				else 
				{
					this_param[l]=false;
					i=parameter_name.length();
				}
			}
		}
		file_.close();
		file_.open(filename_.c_str(), std::ios::out | std::ios::trunc);
		
		
		bool writen;
		for(l=0;l<line_number;l++)
		{
			if(!this_param[l])file_<<lines[l]<<endl;
			else 
			{
				file_<< parameter_name <<'='<<parameter<<endl;
				writen = true;
			}
		}
		
		if(!writen)file_<< parameter_name <<'='<<parameter<<endl;
			
		
		file_.close();
		file_.clear();

    }

#ifndef SERIAL
	parallel.barrier();
#endif
	
}


//SEARCH=====================================
bool SettingsFile::search(const std::string searchString)
{
  unsigned int i=0;
  char c;
  //Set to beginning of file
  stream_.clear(); //clear any errors from having failed a previous search
  stream_.seekg(0);
 
  //Search
  while(stream_.good() && i<searchString.length())
    {
      c=stream_.get();
      if(c==searchString[i]){i++;}
      else{i=0;}
    }
  return stream_.good();
}
}

// Copyright (C) 2017 by SIB Swiss Institute of Bioinformatics, Julien Dorier and Dimos Goundaroulis.
// 
// This file is part of project Knoto-ID.
// 
// Knoto-ID is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
// 
// Knoto-ID is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
// 
// You should have received a copy of the GNU General Public License
// along with Knoto-ID.  If not, see <http://www.gnu.org/licenses/>.

#include "KnottedCore.hh"


using namespace std;

KnottedCore::KnottedCore(
    const std::vector<double>& x_coords,
    const std::vector<double>& y_coords,
    const std::vector<double>& z_coords,
    const std::string& formattedData,
    const std::vector<std::vector<double> > projectionlist_projections,
    const std::vector<double> projectionlist_weights,
    const bool& all_chains,
    const std::string& closure_method,
    const bool& planar,
    const bool& arrow_polynomial,
    const bool& cyclic,
    const bool& cyclic_input,
    const bool& kc_search_path,
    const bool& reduction_3d,
    const bool& simplify_diagram,
    const time_t timeout,
    const std::string& jones_method,
    const long& max_nb_random_moves_III,
    const long& max_nb_unsuccessfull_random_moves_III,
    const bool& debug
):x_coords(x_coords), y_coords(y_coords), z_coords(z_coords), formattedData(formattedData),
   projectionlist_projections(projectionlist_projections), 
  projectionlist_weights(projectionlist_weights), all_chains(all_chains),  
  closure_method(closure_method), planar(planar), arrow_polynomial(arrow_polynomial),
   cyclic(cyclic), cyclic_input(cyclic_input), kc_search_path(kc_search_path), reduction_3d(reduction_3d), 
   simplify_diagram(simplify_diagram), timeout(timeout), jones_method(jones_method), 
   max_nb_random_moves_III(max_nb_random_moves_III), max_nb_unsuccessfull_random_moves_III(max_nb_unsuccessfull_random_moves_III), 
   debug(debug){
}


std::tuple< std::map<std::string, std::vector<std::string>>, std::map<std::string, std::vector<std::string>>, std::map<std::string, std::vector<std::string>>> KnottedCore::Run()
{
    std::istringstream names_db(formattedData);
    map<string,string> map_jones_to_name;
    unsigned long nb_projections=projectionlist_projections.size();

    string line;
    long line_no=0;
    vector<string> name_jones_tmp;
    const boost::regex validvariable("^(A|v|[Lmwpq][0-9]+)$");
      while (names_db.good())
	{
	  line_no++;
	  getline(names_db, line);
	  name_jones_tmp=split_string(line,"\t");
	  if(name_jones_tmp.size()==2)
	    {
	      Polynomial p;
	      if(!p.load_from_string(name_jones_tmp[1]))
		{
		  cerr<<"*********************************************************"<<endl;
		  cerr<<"ERROR names_db:"<<endl;
		  cerr<<"Invalid polynomial (line "<<line_no<<"):"<<endl;
		  cerr<<name_jones_tmp[1]<<endl;
		  cerr<<"*********************************************************"<<endl;
		  exit(1);
		}		    
	      //check variable names
	      vector<string> var_names=p.get_variable_names();
	      for(int i=0;i<var_names.size();i++)
		{
		  if(!regex_match(var_names[i],validvariable,boost::match_default))
		    {
		      cerr<<"*********************************************************"<<endl;
		      cerr<<"Invalid polynomial:"<<endl;
		      cerr<<"Invalid variable name \""<<var_names[i]<<"\". Variable names should be"<<endl;
		      cerr<<" \"A\",\"v\",\"L1\",\"L2\",...,\"m1\",\"m2\",...,\"w1\",\"w2\",..."<<endl;
		      cerr<<" \"p1\",\"p2\",...,\"q1\",\"q2\",..."<<endl;
		      cerr<<"*********************************************************"<<endl;
		      exit(1);		  
		    }
		}
        map<string,string>::iterator it_map;
	      it_map=map_jones_to_name.find(p.to_string());
	      if(it_map== map_jones_to_name.end())
		{
		  map_jones_to_name[p.to_string()]=name_jones_tmp[0];
		}
	      else
		{
		  it_map->second=it_map->second+"|"+name_jones_tmp[0];
		}
	    }
	  else if(name_jones_tmp.size()!=0)
	    {
	      cerr<<"*********************************************************"<<endl;
	      cerr<<"ERROR names_db (line "<<line_no<<"): invalid file format."<<endl;
	      cerr<<line<<endl;
	      cerr<<"*********************************************************"<<endl;
	      exit(1);
	    }
	}
      

      //add special cases
      map_jones_to_name["TIMEOUT"]="TIMEOUT";
      map_jones_to_name["failed_projection"]="failed_projection";


    Polygon polygon(x_coords,y_coords,z_coords, cyclic_input, debug);
    //i1=start index, l=length i2=i1+l=end index
    //moves: i1++,i1--,i2++,i2--
    // <=> i1++,i1--,l++,l--
    //first phase: decrease l until changing type
    //second phase: follow boundary
    int phase=1;
    int N=polygon.get_nb_points();
    int lmax=N-1;
    if(cyclic_input)//for cyclic input curve, l>=N => full CLOSED chain. (trick to avoid being locked)
        lmax=N+1;
    int i1=0;
    int l=lmax;
    int i1_last,l_last;
    int i1_start=i1,l_start=l,direction_start,total_rotation=0;//stop when reaching again this point
    //directions in (i1,l) plane: 0-> (1,-1), 1-> (0,1), 2-> (-1,0), 3-> (0,-1)
    //direction are with respect to (i1_last,l_last)
    int direction=3;//last direction
    int direction_test=3;
    vector<pair<int,int> > directions;
    directions.push_back(make_pair(1,-1));
    directions.push_back(make_pair(0,1));
    directions.push_back(make_pair(-1,1));
    directions.push_back(make_pair(0,-1));
    string jones_initial="";
    map<pair<int,int>,pair<string,double> > pos_to_jones;//pos_to_jones[pair(i1,l)]=pair(jones,frequ);
    //knotted core
    set<pair<int,int> > knotted_cores;
    int best_l=l;
    knotted_cores.insert(make_pair(i1,l));
    int nb_steps=0;
    bool flag_timeout_during_search=false;

    auto mod = [](int x, int N) {
        if(x >= 0)
            return x % N;
        else 
            return (x + N * (2 - x / N)) % N;
    };

    // Output the results
    std::map<std::string, std::vector<std::string>> file_search;

        file_search["index_first"];
        file_search["index_last"];
        file_search["length"];
        file_search["frequency"];

        if(cyclic)
            file_search["knot_type"];
        else
            file_search["knotoid_type"];

        file_search["polynomial"];
  // std::cout<<"\r"<<"Searching for Knotted Core."<<" ";
  while(true)
    {
      // cerr<<"\r"<<" testing subchain "<<mod(i1,N)<<"-"<<mod(i1+l,N)<<"      ";
      string jones;
      double frequ;
      if((l<=0)||(l>lmax)||(cyclic_input==false&&i1<0)||(cyclic_input==false&&i1+l>=N))//i1>=N
	{
	  jones="OUT_OF_BOUNDS";
	  frequ=1;
	}
      else
	{
	  map<pair<int,int>,pair<string,double> >::iterator it_jones;
	  it_jones = pos_to_jones.find(make_pair(mod(i1,N),l));      
	  if (it_jones!= pos_to_jones.end())//jones already evaluated
	    {
	      jones=it_jones->second.first;
	      frequ=it_jones->second.second;
	    }
	  else//evaluate jones
	    {
	      if(cyclic_input&&l>=N)
		{//evaluate jone for full polygon, cyclic
		  jones=get_jones(polygon,frequ,1,"direct");
		}
	      else
		{
		  Polygon polygontmp=polygon.get_polygon(mod(i1,N),mod(i1+l,N),cyclic);
		  jones=get_jones(polygontmp,frequ,nb_projections,closure_method);
		}
	    //   if(jones=="TIMEOUT")
		// {
		//   cerr<<endl;
		//   cerr<<"*********************************************************"<<endl;
		//   cerr<<"WARNING: timeout during evaluation of the polynomial."<<endl;
		//   cerr<<"*********************************************************"<<endl;
		//   flag_timeout_during_search=true;
		// }
	      if(pos_to_jones.size()==0)
      jones_initial=jones;
	      pos_to_jones[make_pair(mod(i1,N),l)]=make_pair(jones,frequ);		  
	    }
	}
      //output
	  string knot_type="UNKNOWN";	        
	      map<string,string>::iterator it=map_jones_to_name.find(jones);
	      if (it != map_jones_to_name.end())
		knot_type=it->second;

        if(jones!="OUT_OF_BOUNDS")
        {
            std::ostringstream oss;

            oss << mod(i1,N);
            file_search["index_first"].push_back(oss.str());
            oss.str(""); // clear the stringstream

            oss << mod(i1+l,N);
            file_search["index_last"].push_back(oss.str());
            oss.str(""); // clear the stringstream

            oss << l+1;
            file_search["length"].push_back(oss.str());
            oss.str(""); // clear the stringstream

            oss << frequ;
            file_search["frequency"].push_back(oss.str());
            oss.str(""); // clear the stringstream

            // if(names_db_filename!="")file_search["knot_type"].push_back(knot_type);

            file_search["polynomial"].push_back(jones);
        }

      //next points
      if(phase==1)
	{
	  if(jones==jones_initial)//valid jones
	    {
	      i1_last=i1;
	      l_last=l;
	      direction=direction_test;
	      direction_test=3;
	    }
	  else//end of phase 1
	    {
	      phase=2;
	      direction=0;
	      direction_test=0;//we already tested direction 3
	      nb_steps=0;
	    }
	}
      else if(phase==2)
	{
	  if(jones==jones_initial)//valid move
	    {
	      nb_steps++;
	      i1_last=i1;
	      l_last=l;
	      if(nb_steps==1)
		{
		  i1_start=i1_last;
		  l_start=l_last;
		  direction_start=direction_test;
		  total_rotation=0;
		}
	      else if(nb_steps>1)
		{
		  if(mod(direction_test-direction,4)==3)//direction_test starts at (direction+3)%4 (i.e. -1) and only increase
		    total_rotation-=1;
		  else
		    total_rotation+=mod(direction_test-direction,4);		    
		}
	      direction=direction_test;
	      direction_test=(direction+3)%4;//last direction
	    }
	  else
	    {
	      direction_test=(direction_test+1)%4;
	    }
	}

      //store best results
      if(jones==jones_initial)
	{
	  if(l<best_l)
	    {
	      best_l=l;
	      knotted_cores.clear();
	      knotted_cores.insert(make_pair(mod(i1,N),l));
	    }
	  else if(l==best_l)
	    {
	      knotted_cores.insert(make_pair(mod(i1,N),l));
	    }
	    
	}
      //move
      i1=i1_last+directions[direction_test].first;
      l=l_last+directions[direction_test].second;
      
      //check for termination
      if(cyclic_input==false&&i1+l>=N)
       	{
       	  break;
      	}

      if(phase==2&&nb_steps>1&&mod(i1_last,N)==i1_start&&l_last==l_start&&direction_start==direction)
       	{
	  if(total_rotation==0)//finished
	    break;
	  else if(total_rotation==-4)//isolated loop.
	    {
	      phase=1;
	      //pick one of the knotted_cores as new starting point
	      i1_last=knotted_cores.begin()->first;
	      l_last=knotted_cores.begin()->second;
	      direction=3;
	      direction_test=3;
	      //move
	      i1=i1_last+directions[direction_test].first;
	      l=l_last+directions[direction_test].second;
	    }
      	}
       if(l>lmax)
       	{
	  cerr<<"OUPS: l>lmax"<<endl;
       	  break;
      	}
    }
  
//   if(output_filename_search!="")
//     if(output_filename_search!="stdout"&&output_filename_search!="-")
//       file_search_tmp.close();

  // cerr<<endl;
  // std::cout<<"Done."<<endl;

// Output knotted core

std::map<std::string, std::vector<std::string>> file_kc;

file_kc["index_first"];
file_kc["index_last"];
file_kc["length"];
file_kc["frequency"];

if(cyclic)
    file_kc["knot_type"];
else
    file_kc["knotoid_type"];
file_kc["polynomial"];

for(set<pair<int,int> >::iterator it=knotted_cores.begin();it!=knotted_cores.end();it++)
{
    int i1=it->first;
    int l=it->second;
    double frequ;
    string jones;
    map<pair<int,int>,pair<string,double> >::iterator it_jones = pos_to_jones.find(make_pair(i1,l));      
    if (it_jones!= pos_to_jones.end())//jones already evaluated
    {
        jones=it_jones->second.first;
        frequ=it_jones->second.second;
    }
    else
    {
        cerr<<"ERROR: should not happen!!!  i1="<<i1<<" l="<<l<<endl;
        exit(1);
    }
    string knot_type="UNKNOWN";	          
    map<string,string>::iterator it_name=map_jones_to_name.find(jones);
    if (it_name != map_jones_to_name.end())
        knot_type=it_name->second;
    if(flag_timeout_during_search)
    {
        knot_type="TIMEOUT:"+knot_type;
        jones="TIMEOUT:"+jones;
    }

    file_kc["index_first"].push_back(std::to_string(i1));
    file_kc["index_last"].push_back(std::to_string(mod(i1+l,N)));
    file_kc["length"].push_back(std::to_string(l+1));
    file_kc["frequency"].push_back(std::to_string(frequ));
    if(cyclic)
        file_kc["knot_type"].push_back(knot_type);
    else
        file_kc["knotoid_type"].push_back(knot_type);
    // file_search["knot_type"].push_back(knot_type);
    file_kc["polynomial"].push_back(jones);
}

    ///////////////////////////
    //Evaluate all subchains
    ///////////////////////////
    if(all_chains==true)
        {
        // std::cout<<"\r"<<"evaluating all subchains"<<" ";     
    std::map<std::string, std::vector<std::string>> file_all;

    file_all["index_first"];
    file_all["index_last"];
    file_all["length"];
    file_all["frequency"];

    if(cyclic)
        file_all["knot_type"];
    else
        file_all["knotoid_type"];

    file_all["polynomial"];

        int N=polygon.get_nb_points();
        for(int i1=0;i1<N;i1++)
        {
        int lmax2;
        if(cyclic_input)
            lmax2=N+1;
        else
            lmax2=N-i1;
        for(int l=1;l<lmax2;l++)
            {
            // cerr<<"\r"<<" subchain "<<i1<<"-"<<mod(i1+l,N)<<"      ";
            Polygon polygontmp=polygon.get_polygon(i1,mod(i1+l,N),cyclic);

            //check if jones already evaluated
            double frequ;
            string jones;
            map<pair<int,int>,pair<string,double> >::iterator it_jones;
            it_jones = pos_to_jones.find(make_pair(i1,l));      
            if (it_jones!= pos_to_jones.end())//jones already evaluated
            {
            jones=it_jones->second.first;
            frequ=it_jones->second.second;
            }
            else//evaluate jones
            {
            if(cyclic_input&&l>=N)
                {//evaluate jone for full polygon, cyclic
                jones=get_jones(polygon,frequ,1,"direct");
                }
            else
                {
                jones=get_jones(polygontmp,frequ,nb_projections,closure_method);
                }		  
            }
            
            //double frequ;
            string knot_type="UNKNOWN";	      
   
            map<string,string>::iterator it=map_jones_to_name.find(jones);
            if (it != map_jones_to_name.end())
                knot_type=it->second;
            file_all["index_first"].push_back(std::to_string(i1));
            file_all["index_last"].push_back(std::to_string(mod(i1+l,N)));
            file_all["length"].push_back(std::to_string(l+1));
            file_all["frequency"].push_back(std::to_string(frequ));
            if(cyclic)
                file_all["knot_type"].push_back(knot_type);
            else
                file_all["knotoid_type"].push_back(knot_type);
            file_all["polynomial"].push_back(jones);
            }
        }

        // cerr<<endl;
        // std::cout<<"Done."<<endl;
        if (kc_search_path==true){
          return std::make_tuple(file_kc, file_all, file_search);
        };
        return std::make_tuple(file_kc, file_all,std::map<std::string, std::vector<std::string>>{});
        };

if (kc_search_path==true){
    return std::make_tuple(file_kc, std::map<std::string, std::vector<std::string>>{}, file_search);
        };
  return std::make_tuple(file_kc, std::map<std::string, std::vector<std::string>>{}, std::map<std::string, std::vector<std::string>>{});
};


string KnottedCore::get_jones(Polygon & polygon,double & frequ, unsigned long nb_projections, std::string closure_method)
{

  double dx,dy,dz,weight;
  bool flag_3d_reduction=true;
  map<string,double> histogram_jones;
  for(unsigned long p=0;p<nb_projections;p++)
    {
      Polygon polygontmp=polygon;
      if(projectionlist_projections.size()>0)
	{
	  dx=projectionlist_projections[p][0];
	  dy=projectionlist_projections[p][1];
	  dz=projectionlist_projections[p][2];
	  weight=projectionlist_weights[p];
	}
      else
	{
	  double u=random01();
	  double v=random01();
	  double theta=2*3.14159265358979*u;
	  double phi=acos(2*v-1);
	  dx=cos(theta)*sin(phi);
	  dy=sin(theta)*sin(phi);
	  dz=cos(phi);
	  weight=1/(double)nb_projections;
	}

      //close polygon
      polygontmp.set_closure(dx,dy,dz,closure_method);
      if(flag_3d_reduction)
	{
	  if(debug)cerr<<"Simplifying 3D curve"<<endl;
	  polygontmp.simplify_polygon(dx,dy,dz);
	  if(debug)cerr<<"3D curve has "<<polygontmp.get_nb_points()<<" vertices"<<endl;
	}

      PlanarDiagram diagram(planar);
      bool flag_valid_projection=true;
      try
	{
	  diagram=polygontmp.get_planar_diagram(dx,dy,dz,planar);
	}
      catch (exception& e)//projection failed  
	{
	  flag_valid_projection=false;
	  cerr<< e.what()<<endl;
	}

      if(flag_valid_projection)
	{
	  diagram.set_debug(debug);
	  if(simplify_diagram)
	    {
	      diagram.simplify();
	      if(max_nb_random_moves_III>0)
		{
		  diagram.simplify_with_random_reidemeister_moves_III(max_nb_random_moves_III,max_nb_unsuccessfull_random_moves_III);
		  diagram.simplify();
		}
	    }
	  //////////evaluate jones//////////////
	  PolynomialInvariant jones(diagram,planar,arrow_polynomial,debug);
	  jones.set_timeout(timeout);
	  Polynomial jones_polynomial;      
	  /////////////jone bgl
	  try
	    {
	      if(jones_method=="simple")
		jones_polynomial=jones.get_polynomial_simple();
	      else if(jones_method=="recursive")
		jones_polynomial=jones.get_polynomial_recursive("default",true);
	      else if(jones_method=="recursive-crossing-order")
		jones_polynomial=jones.get_polynomial_recursive("crossing_order",true);
	      else if(jones_method=="recursive-arc-order")
		jones_polynomial=jones.get_polynomial_recursive("arc_order",true);
	      else if(jones_method=="recursive-region-order")
		jones_polynomial=jones.get_polynomial_recursive("region_order",true);
	      else
		{
		  cerr<<"********************************************************"<<endl;	  
		  cerr<<"ERROR: --polynomial-method="<<jones_method<<" not implemented"<<endl;
		  cerr<<"********************************************************"<<endl;	  
		  exit(1);
		}
	      histogram_jones[jones_polynomial.to_string()]+=weight;
	    }
	  catch (exception& e)//timeout  
	    {
	      histogram_jones["TIMEOUT"]+=weight;
	    }

	}
      else//flag_valid_projection==false
	{
	  histogram_jones["failed_projection"]+=weight;
	}
      
    }

  //find max frequency
  string jones="";
  frequ=0;
  for(map<string,double>::iterator it=histogram_jones.begin();it!=histogram_jones.end();it++)
    {
      if(it->second >= frequ)
	{
	  frequ=it->second;
	  jones=it->first;
	}
    }

  return jones;
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

vector<string> KnottedCore::split_string(string str,string sep_list)
{
  vector<string> result;
  string strtmp;
  size_t pos1=0,pos2,postmp;
  pos2=string::npos;
  for(int i=0;i<sep_list.length();i++){
    postmp=str.find(sep_list[i]);
    if(postmp<pos2)
      pos2=postmp;
  }
   
  while( pos2!=str.npos){
    strtmp=str.substr (pos1,pos2-pos1);
    if(strtmp.size()>0)
      result.push_back(strtmp);
    pos1=pos2+1;
     
    pos2=string::npos;
    for(int i=0;i<sep_list.length();i++){
      postmp=str.find(sep_list[i],pos1);
      if(postmp<pos2)
 	pos2=postmp;
    }
  }
  strtmp=str.substr (pos1,pos2-pos1);
  if(strtmp.size()>0)
    result.push_back(strtmp);
  return result;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////
vector<string> KnottedCore::split_input(string str)
{
  vector<string> result;
  string strtmp;
  string sep="PD[";
  int pos1=0,pos2=str.find(sep);
  
  while( pos2!=str.npos){
    if(pos1>0)//first substring is ignored
      {
	strtmp=str.substr (pos1-1,pos2-pos1+1);
	if(strtmp.size()>0)
	  result.push_back(strtmp);
      }
    pos1=pos2+1;
    pos2=str.find(sep,pos1);    
  }
  if(pos1>0)
    {
      strtmp=str.substr (pos1-1,pos2-pos1+1);
      if(strtmp.size()>0)
	result.push_back(strtmp);
    }
  return result;
}
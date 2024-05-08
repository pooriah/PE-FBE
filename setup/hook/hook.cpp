#include <iostream>
#include <fstream>

#include <kernel/hydro.h>
#include <util/format.h>
#include <util/posthook.h>
#include <ctime>


#include "func/init_bc.h"

// Data to be shared among all local blocks
struct SharedState {
	int group=-1;

	std::string key="electro_neumann";
	
	//bubble_generation
	std::vector<std::vector<double>> gentable,surfacetable; //x0,zy0,z0,lx,ly,lz
	double minsurf,maxsurf;
	int dimsurf;
	
	int total;
	int iterator;
	
	//initial bubbles
	int done=0;
	int n_init_bub=0;
};

template <class T>
T dist(T x1,T y1,T x2,T y2)
{
	return(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)));
}

template <class T>
T dist(T x1,T y1,T z1,T x2,T y2,T z2)
{
	return(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2)));
}

template <class M>
void InitHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  auto& var = hydro->var;
  auto& m = hydro->m;
  auto sem = m.GetSem();
  const Scal bubgen_seq(var.Int["enable_bubgen_seq"]);
  const Scal onsurface_generation(var.Int["onsurface_generation"]);
  if (sem()) {
    if (m.IsLead()) { // Executed only on the lead block, once in each rank.
                      // There is only one lead block in each rank.
      // Allocate memory for the shared state and save pointer in the solver.
      hydro->par_.ptr = new SharedState();
      auto* shared = static_cast<SharedState*>(hydro->par_.ptr);

    }
    if (m.IsRoot()) { // Executed only on the root block, once over all ranks.
                      // There is only one root block over all ranks.
      auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
      
      shared->iterator=0;
      shared->total=-1;
      
      if(bubgen_seq)
      {
	      std::stringstream ss(var.String["bubgen_path"]);
	      std::string temp;
	      while(std::getline(ss,temp,'\n'))
	      {
	      	shared->total++;
	      }
	      
	      size_t pos;
	      std::stringstream sgen(var.String["generation_region"]);
	      std::getline(sgen,temp,'\n');
	      while(std::getline(sgen,temp,'\n'))
	      {
	      	pos=temp.find(" ",0);
	      	shared->gentable.push_back(std::vector<double>());
	      	for(int i=0;i<6;i++)
	      	{
	      		shared->gentable[shared->gentable.size()-1].push_back(std::stod(temp.substr(pos,temp.find(" ",pos+1)-pos-1)));
	      		pos=temp.find(" ",pos+1);
	      	}

	      }
      }
      if(bubgen_seq && onsurface_generation)
      {      
      	      std::string temp;
	      size_t pos;
	      std::stringstream sgen(var.String["generation_electrode"]);
	      std::getline(sgen,temp,'\n');
	      while(std::getline(sgen,temp,'\n'))
	      {
	      	pos=temp.find(" ",0);
	      	shared->surfacetable.push_back(std::vector<double>());
	      	for(int i=0;i<6;i++)
	      	{
	      		shared->surfacetable[shared->surfacetable.size()-1].push_back(std::stod(temp.substr(pos,temp.find(" ",pos+1)-pos-1)));
	      		pos=temp.find(" ",pos+1);
	      	}

	      }
	      
	      int dim=-1;
	      for(int i=3;i<6;i++)
	      {
	      	if(shared->surfacetable[0][i]==0){dim=i-3;break;}
	      }
	      
	      try
	      {
	      	if(dim<0){throw(dim);}
	      }
	      catch(int error)
	      {
	      	std::cout<<"dim is negative: "<<dim<<"\n";
	      }
	      
	      shared->dimsurf=dim;
	      
	      shared->minsurf=shared->surfacetable[0][dim];
	      shared->maxsurf=shared->surfacetable[0][dim];
	      
	      for(int i=1;i<shared->surfacetable.size();i++)
	      {
	      	if(shared->minsurf>shared->surfacetable[i][dim]){shared->minsurf=shared->surfacetable[i][dim];}
	      	if(shared->maxsurf<shared->surfacetable[i][dim]){shared->maxsurf=shared->surfacetable[i][dim];}
	      }
	      	
      }
      	
      
      std::cout << util::Format(
          "InitHook t={:} total={:}\n\n",hydro->fs_->GetTime(),shared->total);
    }
    
    //Initial number of bubbles
    {
    	auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
    	shared->n_init_bub=var.Int["n_init_bub"];
    }
    
    
  }

  if(sem())
  {
      auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
      hydro->me_group_.LoopPairs([&](auto cf_group) {
        auto cf = cf_group.first;
        size_t group = cf_group.second;
        auto& custom = hydro->bc_group_custom_[group];
        auto getptr = [&](std::string key) -> size_t* {
          auto it = custom.find(key);
          if (it != custom.end() && (it->second)>0) {
            return& group;
          }
          return nullptr;
        };

        if(auto* group=getptr(shared->key))
        {
            shared->group=*group;
        }
      });
      if(m.IsRoot())
      {
	if(m.IsRoot())
	{
		std::ofstream file;
		file.open("electro_hist.txt");
		file<<"time\tArea\tCurrent\tPotential\n";
		file.close();
	}
      }
  }
}

template <typename T> std::string type_name();

template <class M>
void StepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto& m = hydro->m;
  auto& eb=*hydro->eb_;
  auto& var = hydro->var;
  auto sem = m.GetSem();
  auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
  const Scal t=hydro->fs_->GetTime();

  // State persistent across all stages (sections `sem()`) of one call
  struct {
  	double current,potential;
  	double area;
  	double count_pot;
  	double test;
  } * ctx(sem);
  auto& vb = *ctx;

  if(sem())
  {
  	if(hydro->dumper_.Try(hydro->st_.t,hydro->st_.dt))
  	{
		auto as=dynamic_cast<typename Hydro<M>::ASVEB*>(hydro->as_.get());
		auto& fcvf=as->GetField();
		auto test=UEmbed<M>::Interpolate(fcvf,hydro->mebc_vfsm_,eb);
		const auto& idc=m.GetIndexCells();
		const auto& idf=m.GetIndexFaces();

	  	auto& cur=hydro->electro_->GetCurrent();
	  	const auto pot=hydro->electro_->GetPotential();
	  	auto mc=hydro->mebc_electro_;

		auto dc_=hydro->me_group_.GetMapCell();
		auto df_=hydro->me_group_.GetMapFace();
		vb.current=0;vb.area=0;vb.count_pot=0;vb.potential=0;
		mc.LoopBCond(eb,[&](auto cf,IdxCell c,auto& bc)
		{
			size_t group=hydro->me_group_[cf];
			auto id=idc.GetMIdx(c);
			auto len=m.GetCellSize();
			if(group==shared->group && m.IsInner(c))
			{
				vb.current+=cur[c].dot(eb.GetNormal(cf))*len[0]*len[1];
				vb.count_pot++;
				vb.potential+=pot[c];
				
				if(fcvf[c]<0.1)
			        	vb.area+=len[0]*len[1];
			}
		});
	    	m.Reduce(&vb.area,Reduction::sum);
		m.Reduce(&vb.current,Reduction::sum);
		m.Reduce(&vb.potential,Reduction::sum);
		m.Reduce(&vb.count_pot,Reduction::sum);
	}
  }

  if(sem())
  {
      if(hydro->dumper_.Try(hydro->st_.t,hydro->st_.dt))
      if(m.IsRoot())
      {
      	//time\tGroup\tArea\tCurrent\tPotential\n
	std::ofstream file;
	file.open("electro_hist.txt",std::ios_base::app);
	file<<t<<"\t"<<vb.area<<"\t"<<vb.current<<"\t"<<vb.potential/vb.count_pot<<"\n";
	std::cout<<t<<"\t"<<vb.area<<"\t"<<vb.current<<"\t"<<vb.potential/vb.count_pot<<"\n";
	file.close();
      }
  }
}

template <class M>
void PreStepHook(Hydro<M>* hydro) {
  using Scal = typename M::Scal;
  using Vect = typename M::Vect;

  auto& m = hydro->m;
  auto& eb=*hydro->eb_;
  auto& var = hydro->var;
  auto sem = m.GetSem();
  auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
  
  const Scal t = hydro->fs_->GetTime();
  auto& stat=hydro->st_;
  const Scal bubgen_seq(var.Int["enable_bubgen_seq"]);
  const Scal onsurface_generation(var.Int["onsurface_generation"]);
  std::stringstream ss(var.String["bubgen_path"]);
  std::string temp,bubble;
  
  using EB = Embed<M>;

  // State persistent across all stages (sections `sem()`) of one call
  struct {
  	//For current
  	double area;
  	
  	//For trajectory
	std::vector<std::string> column_names;
	std::vector<Scal> row_colors;
	std::vector<std::vector<Scal>> table;

	//For bubble generation
	std::shared_ptr<FieldCell<Scal>> fcvf; // volume fraction
	Vars var;
	bool verbose = false;
  } * ctx(sem);
  auto& vb = *ctx;
  
  if(sem.Nested("trajectory"))
  {
	auto& k = *hydro;
	auto plic=hydro->as_->GetPlic();
	if(hydro->eb_)
	{
		CalcTraj<EB>(
		  *k.eb_, k.layers, plic.vfcu, plic.vfccl_stat, plic.vfcim,
		  k.fs_->GetPressure(), k.fs_->GetVelocity(), vb.column_names,
		  vb.row_colors, vb.table);
	}
	else
	{
		CalcTraj<M>(
		  m, k.layers, plic.vfcu, plic.vfccl_stat, plic.vfcim,
		  k.fs_->GetPressure(), k.fs_->GetVelocity(), vb.column_names,
		  vb.row_colors, vb.table);
	}
  }
  
  if(sem())
  {
  }
  	
  
  if(bubgen_seq)
  {
  
  auto& par = *ctx;
  
  if (hydro->bubgen_.Try(stat.t, stat.dt)) {
    if (sem("as-bubgen-var"))
    {
	if(m.IsRoot())
	{
	  	std::cout<<"Iterator is: "<<shared->iterator<<"\n";
		std::getline(ss,bubble);
		for(int i=0;i<=shared->iterator;i++)
		{
			std::getline(ss,temp);
		}
		auto as=dynamic_cast<typename Hydro<M>::ASVEB*>(hydro->as_.get());
	auto& fcvf=as->GetField();
	const auto& idc=m.GetIndexCells();
	const auto& idf=m.GetIndexFaces();

  	auto& cur=hydro->electro_->GetCurrent();
  	auto mc=hydro->mebc_electro_;

	vb.area=0;
	mc.LoopBCond(eb,[&](auto cf,IdxCell c,auto& bc)
	{
		size_t group=hydro->me_group_[cf];
		auto id=idc.GetMIdx(c);
        	auto len=m.GetCellSize();
		if(group==shared->group && m.IsInner(c))
		{
			if(fcvf[c]<0.1)
                		vb.area+=len[0]*len[1];
		}
	});
		const Scal lx(var.Double["extent"]),ly=lx;
		std::srand(time(nullptr));
		
		Scal off=1.4;
		std::string out="",out_temp;
		size_t pos;
		
		std::vector<size_t> ind;
		Scal x0,y0,z0,x,y,z;
		int nx,ny,nz;
		std::string sub=" ",shape="sphere";
		double rr;
		Scal temp_dst,dx;
		size_t nb;
		Scal rmin(var.Double["rmin"]);
		
		for(int cc=0;cc<shared->gentable.size();cc++)
		{
			if(shared->done==0 && shared->n_init_bub>0)
			{
				nb=shared->n_init_bub;
				
				if(cc==(shared->gentable.size()-1))
				{
					shared->n_init_bub=0;
				}
				
				std::cout<<"Init bub "<<nb<<"\n";

			}
			else
			{
				nb=1;
				std::cout<<"---- Here\n";
			}
			
			rr=std::stod(temp.substr(temp.find_last_of(sub)+1))/pow(pow(var.Double["gen_ratio"],cc),1.0/3.0);	
			ind.clear();
			while(rr>rmin)
			{
				std::cout<<"####\t"<<nb<<"\t"<<rr<<"\n";
				
				dx=rr/2.0;
				if(shared->done==0)
				{
					dx=(rr*3.0);
					if(cc==(shared->gentable.size()-1))
					{
						shared->done=1;
					}
				}
				
				nx=(shared->gentable[cc][3]-2*rr*off)/dx+1;
				ny=(shared->gentable[cc][4]-2*rr*off)/dx+1;
				nz=(shared->gentable[cc][5]-2*rr*off)/dx+1;
				x0=(shared->gentable[cc][3]-(nx-1)*dx)/2+shared->gentable[cc][0];
				y0=(shared->gentable[cc][4]-(ny-1)*dx)/2+shared->gentable[cc][1];
				z0=(shared->gentable[cc][5]-(nz-1)*dx)/2+shared->gentable[cc][2];
				
				if(onsurface_generation){
				if(shared->surfacetable[cc][3]==0)
				{
					nx=1;
					if((abs(shared->surfacetable[cc][0]-shared->minsurf)/(shared->surfacetable[cc][0]+1e-3))<1e-3)
					{
						x0=shared->surfacetable[cc][0]+rr*off;
					}
					else
					{
						x0=shared->surfacetable[cc][0]-rr*off;
					}
				}
				else if(shared->surfacetable[cc][4]==0)
				{
					ny=1;
					std::cout<<shared->surfacetable[cc][1]<<"\t"<<shared->minsurf<<" ----\n";
					if((abs(shared->surfacetable[cc][1]-shared->minsurf)/(shared->surfacetable[cc][1]+1e-3))<1e-3)
					{
						y0=shared->surfacetable[cc][1]+rr*off;
						std::cout<<"1-"<<y0<<"\n";
					}
					else
					{
						y0=shared->surfacetable[cc][1]-rr*off;
						std::cout<<"2-"<<y0<<"\n";
					}
				}
				else if(shared->surfacetable[cc][5]==0)
				{
					nz=1;
					if((abs(shared->surfacetable[cc][2]-shared->minsurf)/(shared->surfacetable[cc][2]+1e-3))<1e-3)
					{
						z0=shared->surfacetable[cc][2]+rr*off;
					}
					else
					{
						z0=shared->surfacetable[cc][2]-rr*off;
					}
				}}
				
	
				
				for(int k=0;k<nz;k++)
				{
					for(int i=0;i<nx;i++)
					{
						for(int j=0;j<ny;j++)
						{
							x=x0+i*dx;
							y=y0+j*dx;
							z=z0+k*dx;
							for(int ii=0;ii<vb.table.size();ii++)
							{
								temp_dst=dist(x,y,z,vb.table[ii][2],vb.table[ii][3],vb.table[ii][4])-rr*off-vb.table[ii][1];
								if(temp_dst<0)
								{
									ind.push_back(k*nx*ny+i*ny+j);
									break;
								}
							}
						}
					}
				}
				
				
				if((nx*ny*nz)<(nb+ind.size()) || nx<=0 || ny<=0 || nz<=0)
				{
					rr=rr/pow(2,1.0/3);
					nb=nb*2;
					ind.clear();
					continue;
				}
				else
				{
					std::vector<int> place;
					size_t place_temp;
					int ind_size;
					
					
					for(int i=0;i<nb;i++)
					{
						while(ind.size()<(nx*ny*nz))
						{	
							
							place_temp=rand()%(nx*ny*nz-ind.size());
							for(int j=0;j<ind.size();j++)
							{
								if(place_temp>=ind[j])
									place_temp++;
							}
							z=(place_temp/(nx*ny))*dx+z0;
							x=((place_temp-nx*ny*(place_temp/(nx*ny)))/ny)*dx+x0;
							y=((place_temp-nx*ny*(place_temp/(nx*ny)))%ny)*dx+y0;
							
							
							
							for(int ii=0;ii<vb.table.size();ii++)
							{
								temp_dst=dist(x,y,z,vb.table[ii][2],vb.table[ii][3],vb.table[ii][4])-rr*off-vb.table[ii][1];
								
								
							}
							
							auto it=find(ind.begin(),ind.end(),place_temp);
							
							
							
							
							ind_size=ind.size();
							for(int k=0;k<place.size();k++)
							{
								if((dist(x,y,z,((place[k]-nx*ny*(place_temp/(nx*ny)))/ny)*dx+x0,((place[k]-nx*ny*(place_temp/(nx*ny)))%ny)*dx+y0,(place[k]/(nx*ny))*dx+z0)-rr*off-rr)<0)
								{
									ind.push_back(place_temp);
									std::sort(ind.begin(),ind.end());
									break;
								}
							}
							if(ind_size==ind.size())
							{
								place.push_back(place_temp);
								break;
							}
						}
					}
					if(place.size()<nb)
					{
						std::cout<<"Not this time\n";
						rr=rr/pow(2,1/3);
						nb=nb*2;
						ind.clear();
						place.clear();
						continue;
					}
					

	
					for(int i=0;i<nb;i++)
					{
						out_temp=temp;
						pos=out_temp.find(shape,0);
						//x_pos
						pos=out_temp.find(sub,pos+1);
						out_temp.replace(pos+1,out_temp.find(sub,pos+1)-pos-1,std::to_string(((place[i]-nx*ny*(place[i]/(nx*ny)))/ny)*dx+x0));
						
						//y_pos
						pos=out_temp.find(sub,pos+1);
						out_temp.replace(pos+1,out_temp.find(sub,pos+1)-pos-1,std::to_string(((place[i]-nx*ny*(place[i]/(nx*ny)))%ny)*dx+y0));
						
						//z_pos
						pos=out_temp.find(sub,pos+1);
						out_temp.replace(pos+1,out_temp.find(sub,pos+1)-pos-1,std::to_string((place[i]/(nx*ny))*dx+z0));
						
						//r
						pos=out_temp.find(sub,pos+1);
						out_temp.replace(pos+1,out_temp.find(sub,pos+1)-pos-1,std::to_string(rr));
						
						out=out+out_temp+"\n";
						
						
					}
					
					break;
				}
			}
				
				
			
			if(rr<=rmin)
			{
				out="";
				break;
			}
		
			
		}
		temp=out;
	
		std::cout<<"++++++++++++++++++++++\n";
		bubble=bubble+"\n"+temp;
		std::cout<<bubble<<"\n";
		
		shared->iterator=(shared->iterator+1)%(shared->total);
		
		std::ofstream file;
		file.open("bub_generation_hist.txt",std::ios_base::app);
		file<<bubble;
		file.close();
	}
	
	par.var.String.Set("init_vf", "list");
	par.var.String.Set("list_path",bubble);
	par.var.Int.Set("dim", var.Int["dim"]);
	par.var.Int.Set("list_ls", var.Int["list_ls"]);
	par.fcvf = std::make_shared<FieldCell<Scal>>(m);
	par.verbose = var.Int("verbose_bubgen", 0);
    }
    if (sem.Nested("as-bubgen-initvf")) {
      InitVf(*par.fcvf, par.var, m, par.verbose);
    }
    if (sem("as-bubgen-apply")) {
      const Scal clnew = hydro->fs_->GetTime();
      auto modify = [fcvf = par.fcvf, clnew](auto& u, auto& cl, auto& eb) {
        for (auto c : eb.AllCells()) {
          if ((*fcvf)[c] > 0) {
            const Scal v = std::min((*fcvf)[c], eb.GetVolumeFraction(c));
            u[c] = std::max(u[c], v);
            cl[c] = clnew;
          }
        }
      };
      auto modifym = [fcvf = par.fcvf, clnew](auto& u, auto& cl, auto, auto& eb) {
        for (auto c : eb.AllCells()) {
          if ((*fcvf)[c] > 0) {
            const Scal v = std::min((*fcvf)[c], eb.GetVolumeFraction(c));
            (*u[0])[c] = std::max((*u[0])[c], v);
            (*cl[0])[c] = clnew;
          }
        }
      };
      if (auto* as = dynamic_cast<typename Hydro<M>::ASVEB*>(hydro->as_.get())) {
        as->AddModifier(modify);
      }
      if (auto* as = dynamic_cast<typename Hydro<M>::ASVMEB*>(hydro->as_.get())) {
        as->AddModifier(modifym);
      }
      if (auto* as = dynamic_cast<typename Hydro<M>::ASV*>(hydro->as_.get())) {
        as->AddModifier(modify);
      }
      if (auto* as = dynamic_cast<typename Hydro<M>::ASVM*>(hydro->as_.get())) {
        as->AddModifier(modifym);
      }
    }
  }
  }

  //find the length or area that is not in contact with bubble
  if(sem())
  {   
	auto as=dynamic_cast<typename Hydro<M>::ASVEB*>(hydro->as_.get());
	auto& fcvf=as->GetField();
	const auto& idc=m.GetIndexCells();
	const auto& idf=m.GetIndexFaces();

  	auto& cur=hydro->electro_->GetCurrent();
  	auto mc=hydro->mebc_electro_;

	vb.area=0;
	mc.LoopBCond(eb,[&](auto cf,IdxCell c,auto& bc)
	{
		size_t group=hydro->me_group_[cf];
		auto id=idc.GetMIdx(c);
        	auto len=m.GetCellSize();
		if(group==shared->group && m.IsInner(c))
		{
			if(fcvf[c]<0.1)
                		vb.area+=len[0]*len[1];
		}
	});

    	m.Reduce(&vb.area,Reduction::sum);
  }

  //modify boundary condition to keep current constant
  if(sem())
  {
    const Scal current(var.Double["current"]);
    double bc_new=current/vb.area*var.Double["resist1"];
    auto mc=hydro->mebc_electro_;

    hydro->me_group_.LoopPairs([&](auto cf_group) {
        auto cf = cf_group.first;
        size_t group = cf_group.second;
        auto nci = hydro->mebc_fluid_[cf].nci;


        auto& custom = hydro->bc_group_custom_[group];
        auto getptr = [&](std::string key) -> Scal* {
          auto it = custom.find(key);
          if (it != custom.end()) {
              if(shared->group==hydro->me_group_[cf])
              {
                  it->second=bc_new;
              }
            return &it->second;
          }
          return nullptr;
        };

        if (auto* dirichlet = getptr("electro_dirichlet")) {
          hydro->mebc_electro_[cf] =
              BCond<Scal>(BCondType::dirichlet, nci, *dirichlet);
        } else if (auto* neumann = getptr("electro_neumann")) {
          hydro->mebc_electro_[cf] = BCond<Scal>(BCondType::neumann, nci, *neumann);
        } else {
          fassert(
              false,
              "Unknown electro conditions for group" + std::to_string(group));
        }
      });
  }
  
  
}

template <class M>
void FinalHook(Hydro<M>* hydro) {
  auto& m = hydro->m;
  auto sem = m.GetSem();
  auto* shared = static_cast<SharedState*>(hydro->par_.ptr);
  if (sem()) {
    if (m.IsRoot()) {
      std::cout << util::Format(
          "FinalHook t={:} foo={:}\n\n", hydro->fs_->GetTime(), shared->iterator);
    }
    if (m.IsLead()) {
      delete shared;
    }
  }
}

using M = MeshCartesian<double, 3>;
template void InitHook(Hydro<M>*);
template void StepHook(Hydro<M>*);
template void PreStepHook(Hydro<M>*);
template void FinalHook(Hydro<M>*);

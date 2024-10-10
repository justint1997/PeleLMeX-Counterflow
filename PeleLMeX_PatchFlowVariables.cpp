
#include "PeleLMeX.H"
using namespace amrex;

void
patchFlowVariables(
  const amrex::Geometry& geom, ProbParm const& lprobparm, amrex::MultiFab& a_mf)
  //const amrex::Geometry& geom, ProbParm const& prob_parm, amrex::MultiFab& a_mf)
  
  
					//int i, int j, int k,
                    //int is_incompressible,
                    //amrex::Array4<amrex::Real> const& state,
                    //amrex::Array4<amrex::Real> const& aux,
                    //amrex::GeometryData const& geomdata,
                    //ProbParm const& prob_parm,
                    //pele::physics::PMF::PmfData::DataContainer const * pmf_data)
{

  amrex::Print() << "\nPatching flow variables..";
  //const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo = geom.ProbLoArray();
  //const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx = geom.CellSizeArray();
  
  //const amrex::Real* prob_lo = geomdata.ProbLo();
  //const amrex::Real* prob_hi = geomdata.ProbHi();
  //const amrex::Real* dx      = geomdata.CellSize();
  
  //const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_lo = geom.ProbLo();
  //const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> prob_hi = geom.ProbHi();
  //const amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> dx      = geom.CellSize();
  
  const amrex::Real* prob_lo = geom.ProbLo();
  const amrex::Real* prob_hi = geom.ProbHi();
  const amrex::Real* dx      = geom.CellSize();

  for (amrex::MFIter mfi(a_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi) {
    const amrex::Box& bx = mfi.tilebox();
    
	auto const& rho_arr = a_mf.array(mfi, DENSITY);
    auto const& rhoY_arr = a_mf.array(mfi, FIRSTSPEC);
    auto const& rhoH_arr = a_mf.array(mfi, RHOH);
	auto const& temp_arr = a_mf.array(mfi, TEMP);
	Real massfrac[NUM_SPECIES] = {0.0}; //In PeleLMeX_Plot.cpp
	//amrex::Real massfrac[NUM_SPECIES] = {0.0}; //In pelelmex_prob.H
    
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) noexcept {
      
	  auto eos = pele::physics::PhysicsType::eos();
      amrex::Real massfrac[NUM_SPECIES] = {0.0};
      amrex::Real sumYs = 0.0;
	  
	  for (int n = 0; n < NUM_SPECIES; n++) {
        massfrac[n] = rhoY_arr(i, j, k, n);
#ifdef N2_ID
        if (n != N2_ID) {
          sumYs += massfrac[n];
        }
#endif
      }
#ifdef N2_ID
      massfrac[N2_ID] = 1.0 - sumYs;
#endif
	  
	  //amrex::Real x[3] = {
        //prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
        //prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
        //prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};
       
	  AMREX_D_TERM(const amrex::Real x = prob_lo[0] + (i+0.5)*dx[0];,
               const amrex::Real y = prob_lo[1] + (j+0.5)*dx[1];,
               const amrex::Real z = prob_lo[2] + (k+0.5)*dx[2];);

      AMREX_D_TERM(const amrex::Real Lx = prob_hi[0] - prob_lo[0];,
               const amrex::Real Ly = prob_hi[1] - prob_lo[1];,
               const amrex::Real Lz = prob_hi[2] - prob_lo[2]);

      AMREX_D_TERM(const amrex::Real xc = prob_lo[0] + Lx/2.0;,
               const amrex::Real yc = prob_lo[1] + Lx/2.0;,
               const amrex::Real zc = prob_lo[2] + Lz/2.0);
	  
	  //amrex::Real x[3] = {
        //prob_lo[0] + static_cast<amrex::Real>(i + 0.5) * dx[0],
        //prob_lo[1] + static_cast<amrex::Real>(j + 0.5) * dx[1],
        //prob_lo[2] + static_cast<amrex::Real>(k + 0.5) * dx[2]};
		
	  //amrex::Real Lx[3] = {
        //prob_hi[0] - prob_lo[0];,
        //prob_hi[1] - prob_lo[1];,
        //prob_hi[2] - prob_lo[2]);

      //amrex::Real xc[3] = {
        //prob_lo[0] + Lx[0]/2.0;,
        //prob_lo[1] + Lx[1]/2.0;,
        //prob_lo[2] + Lz[2]/2.0);
	  
      //amrex::ignore_unused(x);
      /*User can define how to patch flow variables here.*/
      //temp_arr(i, j, k) = lprobparm.T_mean;
	  
	  amrex::Real radiusSq = AMREX_D_TERM(  (x-xc) * (x-xc),
                                           + (y-yc) * (y-yc),
                                           + (z-zc) * (z-zc));
	  //amrex::Real radiusSq = AMREX_D_TERM(  (x[0]-xc[0]) * (x[0]-xc[0]),
                                           //+ (y[1]-yc[1]) * (y[1]-yc[1]),
                                           //+ (z[2]-zc[2]) * (z[2]-zc[2]));									   
      amrex::Real radius = std::sqrt(radiusSq);
      //amrex::Real mixingWidth = 0.1*prob_parm.ignitSphereRad;
	  amrex::Real mixingWidth = 0.1*lprobparm.ignitSphereRad;
      //amrex::Real mixingFunction = 0.5 * ( 1.0 + std::tanh((prob_parm.ignitSphereRad - radius)/mixingWidth));
	  amrex::Real mixingFunction = 0.5 * ( 1.0 + std::tanh((lprobparm.ignitSphereRad - radius)/mixingWidth));
	  //state(i,j,k,TEMP) = mixingFunction * prob_parm.ignitT + (1.0 - mixingFunction) * prob_parm.T_inert;
	  temp_arr(i, j, k) = mixingFunction * lprobparm.ignitT + (1.0 - mixingFunction) * lprobparm.T_inert;
	  
	  if (radius <= lprobparm.ignitSphereRad) {
	  
	    massfrac[NC12H26_ID] = 0;
        massfrac[H_ID] = 0.000015;
        massfrac[O_ID] = 0.000224;
        massfrac[OH_ID] = 0.002595;
        massfrac[HO2_ID] = 0.000002;
        massfrac[H2_ID] = 0.000201;
        massfrac[H2O_ID] = 0.08308; //Use equilibrium
        massfrac[H2O2_ID] = 0;
        massfrac[O2_ID] = 0.008535; //Use equilibrium
        massfrac[CH3_ID] = 0;
        massfrac[CH4_ID] = 0;
        massfrac[CH2O_ID] = 0;
        massfrac[CO_ID] = 0.01447; //Use equilibrium
        massfrac[CO2_ID] = 0.17207; //Use equilibrium
        massfrac[C2H2_ID] = 0;
        massfrac[C2H4_ID] = 0;
        massfrac[C2H6_ID] = 0;
        massfrac[CH2CHO_ID] = 0;
        massfrac[aC3H5_ID] = 0;
        massfrac[C3H6_ID] = 0;
        massfrac[C2H3CHO_ID] = 0;
        massfrac[C4H7_ID] = 0;
        massfrac[C4H81_ID] = 0;
        massfrac[C5H9_ID] = 0;
        massfrac[C5H10_ID] = 0;
        massfrac[C6H12_ID] = 0;
        massfrac[C7H14_ID] = 0;
        massfrac[C8H16_ID] = 0;
        massfrac[C9H18_ID] = 0;
        massfrac[PXC9H19_ID] = 0;
        massfrac[C10H20_ID] = 0;
        massfrac[C12H24_ID] = 0;
        massfrac[C12H25O2_ID] = 0;
        massfrac[OC12H23OOH_ID] = 0;
        massfrac[N2_ID] = 0.718808; //Use equilibrium
        }
	  
	  
    });
  }

  amrex::Print() << "Done\n";
}

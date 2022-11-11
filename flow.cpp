//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file flow.cpp
//! \brief Shearing wave problem generator with added routines
//======================================================================================

// C headers

// C++ headers
#include <cmath>      // sqrt()
#include <fstream>    // ofstream
#include <iomanip>    // setprecision
#include <iostream>   // cout, endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()
#include <sstream>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../orbital_advection/orbital_advection.hpp"
#include "../outputs/outputs.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

#if MAGNETIC_FIELDS_ENABLED
#error "This problem generator requires does not work with MHD."
#endif

namespace {
Real iso_cs, gm1, d0, p0, mTrack;
Real amp; // amplitude
int nwx, nwy; // Wavenumbers
Real x1size,x2size,x3size;
int ipert; // initial pattern
Real qshear, Omega0;
Real hst_dt, hst_next_time;
bool error_output;

Real Historydvyc(MeshBlock *pmb, int iout);
Real Historyvxs(MeshBlock *pmb, int iout);
Real Historydvys(MeshBlock *pmb, int iout);
} // namespace

void sink(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

void suck(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

void fill(MeshBlock *pmb, const Real time, const Real dt,
	     const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
	     const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
	     AthenaArray<Real> &cons_scalar);

void lowerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
             FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju,
	     int kl, int ku, int ngh);

void upperX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
	     FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju,
	     int kl, int ku, int ngh);

/*
void bbCool(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);
*/
Real CustomTimeStep(MeshBlock *pmb);

//======================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief Init the Mesh properties
//======================================================================================

void Mesh::InitUserMeshData(ParameterInput *pin) {
  // custom add-ons
  // specify routine inclusion
  if (pin->GetReal("custom","sinkOn") == 0) {
    EnrollUserExplicitSourceFunction(sink);
  }
  if (pin->GetReal("custom","suckOn") == 0) {
    EnrollUserExplicitSourceFunction(suck);
  }
  if (pin->GetReal("custom","fillLength") > 0) {
    EnrollUserExplicitSourceFunction(fill);
  }
  if (pin->GetReal("custom","fixdt") > 0) {
    EnrollUserTimeStepFunction(CustomTimeStep);
  }

  // specify data fields
  AllocateRealUserMeshDataField(1);
  ruser_mesh_data[0].NewAthenaArray(11);
  ruser_mesh_data[0](0) = pin->GetReal("custom","G");
  ruser_mesh_data[0](1) = pin->GetReal("custom","m"); // assignment of m by real works here, in Mesh function
  ruser_mesh_data[0](2) = pin->GetReal("custom","rs");
  ruser_mesh_data[0](3) = pin->GetReal("custom","fixdt");
  ruser_mesh_data[0](4) = pin->GetReal("custom","consMass");
  ruser_mesh_data[0](5) = pin->GetReal("custom","fillLength");
  ruser_mesh_data[0](6) = pin->GetReal("custom","rho0");
  ruser_mesh_data[0](7) = pin->GetReal("orbital_advection","qshear") * pin->GetReal("orbital_advection","Omega0");
  ruser_mesh_data[0](8) = pin->GetReal("hydro","iso_sound_speed");
  ruser_mesh_data[0](9) = pin->GetReal("hydro","gamma");
  ruser_mesh_data[0](10) = pin->GetReal("custom","gmmSink");

  // enroll boundary conditions
  //EnrollUserBoundaryFunction(BoundaryFace::inner_x2, lowerX2);
  //EnrollUserBoundaryFunction(BoundaryFace::outer_x2, upperX2);

  if (!shear_periodic) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
        << "This problem generator requires shearing box."   << std::endl;
    ATHENA_ERROR(msg);
  }

  if (mesh_size.nx2 == 1) {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
        << "This problem does NOT work on a 1D grid." << std::endl;
    ATHENA_ERROR(msg);
  }

  // read ipert parameter
  //d0 = 1.0;
  d0 = pin->GetReal("custom","rho0");
  p0 = 1e-6;
  if (NON_BAROTROPIC_EOS) {
    gm1 = (pin->GetReal("hydro","gamma") - 1.0);
    iso_cs = std::sqrt((gm1+1.0)*p0/d0);
  } else {
    iso_cs = pin->GetReal("hydro","iso_sound_speed");
    p0 = d0*SQR(iso_cs);
  }
  ipert = pin->GetInteger("problem","ipert");
  x1size = mesh_size.x1max - mesh_size.x1min;
  x2size = mesh_size.x2max - mesh_size.x2min;
  x3size = mesh_size.x3max - mesh_size.x3min;

  // shearing box parameters
  qshear = pin->GetReal("orbital_advection","qshear");
  Omega0 = pin->GetReal("orbital_advection","Omega0");

  if (ipert == 3) {
    amp = pin->GetReal("problem","amp");
    nwx = pin->GetInteger("problem","nwx");
    nwy = pin->GetInteger("problem","nwy");
    if (nwx == 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
          << "Parameter nwx must be non-zero." << std::endl;
      ATHENA_ERROR(msg);
    }
    error_output = pin->GetOrAddBoolean("problem","error_output",false);
    if (error_output) {
      // allocateDataField
      AllocateRealUserMeshDataField(2);
      ruser_mesh_data[0].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);
      ruser_mesh_data[1].NewAthenaArray(mesh_size.nx3, mesh_size.nx2, mesh_size.nx1);

      // read history output timing
      InputBlock *pib = pin->pfirst_block;
      while (pib != nullptr) {
        if (pib->block_name.compare(0, 6, "output") == 0) {
          OutputParameters op;
          std::string outn = pib->block_name.substr(6);
          op.block_number = atoi(outn.c_str());
          op.block_name.assign(pib->block_name);
          op.next_time = pin->GetOrAddReal(op.block_name,"next_time", time);
          op.dt = pin->GetReal(op.block_name,"dt");
          op.file_type = pin->GetString(op.block_name,"file_type");
          if (op.file_type.compare("hst") == 0) {
            hst_dt = op.dt;
            hst_next_time = op.next_time;
          }
        }
        pib = pib->pnext;
      }

      // allocate User-defined History Output
      AllocateUserHistoryOutput(3);
      EnrollUserHistoryOutput(0, Historydvyc, "dvyc",
                              UserHistoryOperation::sum);
      EnrollUserHistoryOutput(1, Historyvxs,  "vxs",
                              UserHistoryOperation::sum);
      EnrollUserHistoryOutput(2, Historydvys, "dvys",
                              UserHistoryOperation::sum);
    }
  } else if (ipert == 1 || ipert == 2) {
    amp = 0.0;
    nwx = 0;
    nwy = 0;
  } else {
    std::stringstream msg;
    msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
        << "This problem requires that ipert is from 1 to 3." << std::endl;
    ATHENA_ERROR(msg);
  }
  return;
}

//======================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  int shboxcoord = porb->shboxcoord;
  int il = is - NGHOST; int iu = ie + NGHOST;
  int jl = js - NGHOST; int ju = je + NGHOST;
  int kl = ks;          int ku = ke;
  if (block_size.nx3 > 1) {
    kl = ks - NGHOST;
    ku = ke + NGHOST;
  }

  if (gid == 0) {
    std::cout << "iso_cs = " << iso_cs << std::endl;
    std::cout << "d0 = " << d0 << std::endl;
    std::cout << "p0 = " << p0 << std::endl;
    std::cout << "ipert  = " << ipert  << std::endl;
    std::cout << "[ssheet.cpp]: [Lx,Ly,Lz] = [" <<x1size <<","<<x2size
              <<","<<x3size<<"]"<<std::endl;
  }

  // calculate wave number just for ipert = 3
  Real kx(0.0), ky(0.0);
  if (ipert==3) {
    kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
    ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  }

  Real x1, x2, rd, rp, rvx, rvy;
  // update the physical variables as initial conditions
  for (int k=kl; k<=ku; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        x1 = pcoord->x1v(i);
        x2 = pcoord->x2v(j);
        rd = d0;
        rp = p0;
        if (ipert == 1) {
          // 1) pure shear bg flow:
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = 0.0;
          phydro->u(IM2,k,j,i) = 0.0;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
          phydro->u(IM3,k,j,i) = 0.0;
        } else if (ipert == 2) {
          // 2) epicyclic oscillation
          if (shboxcoord == 1) { // x-y shear
            rvx = 0.1*iso_cs;
            rvy = 0.0;
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = -rd*rvy;
            if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
            phydro->u(IM3,k,j,i) = 0.0;
          } else { // x-z plane
            rvx = 0.1*iso_cs;
            rvy = 0.0;
            phydro->u(IDN,k,j,i) = rd;
            phydro->u(IM1,k,j,i) = rd*rvx;
            phydro->u(IM2,k,j,i) = 0.0;
            phydro->u(IM3,k,j,i) = -rd*(rvy+qshear*Omega0*x1);
          }
        } else if (ipert == 3) {
          // 3) JG HD shwave test
          rvx = amp*iso_cs*std::cos(kx*x1 + ky*x2);
          rvy = amp*iso_cs*(ky/kx)*std::cos(kx*x1 + ky*x2);
          phydro->u(IDN,k,j,i) = rd;
          phydro->u(IM1,k,j,i) = -rd*rvx;
          phydro->u(IM2,k,j,i) = -rd*rvy;
          if(!porb->orbital_advection_defined)
            phydro->u(IM2,k,j,i) -= rd*qshear*Omega0*x1;
          phydro->u(IM3,k,j,i) = 0.0;
        } else {
          std::stringstream msg;
          msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator" << std::endl
              << "Shearing wave sheet ipert=" << ipert << " is unrecognized" << std::endl;
          ATHENA_ERROR(msg);
        }
        if (NON_BAROTROPIC_EOS) {
          phydro->u(IEN,k,j,i) = rp/gm1 + 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                               SQR(phydro->u(IM2,k,j,i)) +
                                               SQR(phydro->u(IM3,k,j,i))
                                              ) / phydro->u(IDN,k,j,i);
        }
      }
    }
  }

  return;
}

//========================================================================================
//! \fn void Mesh::UserWorkInLoop()
//  \brief Function called once every time step for user-defined work.
//========================================================================================

void Mesh::UserWorkInLoop() {
  bool flag = false;
  // check output
  Real present_time = time + dt;
  int cncycle = ncycle + 1;
  if (error_output) {
    if ((present_time < tlim) && (nlim < 0 || cncycle < nlim)
        && (present_time > hst_next_time)) {
      flag = true;
      hst_next_time += hst_dt;
    }
    if ((present_time >= tlim) || (nlim >= 0 && cncycle >= nlim)) {
      flag = true;
    }
  }
  // calculate vxs, dvys
  if (flag) {
    AthenaArray<Real> &vxs = ruser_mesh_data[0];
    AthenaArray<Real> &dvys = ruser_mesh_data[1];
    Real qom = qshear*Omega0;
    Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
    Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
    kx += qom*present_time*ky;

    // initialize vs
    for (int k=0; k<mesh_size.nx3; k++) {
      for (int j=0; j<mesh_size.nx2; j++) {
        for (int i=0; i<mesh_size.nx1; i++) {
          vxs(k,j,i) = 0.0;
          dvys(k,j,i) = 0.0;
        }
      }
    }
    for (int bn=0; bn<nblocal; ++bn) {
      MeshBlock *pmb = my_blocks(bn);
      LogicalLocation &loc = pmb->loc;
      if (loc.level == root_level) { // root level
        for (int k=pmb->ks; k<=pmb->ke; k++) {
          for (int j=pmb->js; j<=pmb->je; j++) {
            for (int i=pmb->is; i<=pmb->ie; i++) {
              int ti = static_cast<int>(loc.lx1)*pmb->block_size.nx1+(i-pmb->is);
              int tj = static_cast<int>(loc.lx2)*pmb->block_size.nx2+(j-pmb->js);
              int tk = static_cast<int>(loc.lx3)*pmb->block_size.nx3+(k-pmb->ks);
              Real x1 = pmb->pcoord->x1v(i);
              Real x2 = pmb->pcoord->x2v(j);
              Real vx = pmb->phydro->w(IVX,k,j,i);
              Real dvy = pmb->phydro->w(IVY,k,j,i);
              if(!pmb->porb->orbital_advection_defined)
                dvy += qom*x1;
              Real vol = pmb->pcoord->GetCellVolume(k,j,i);
              Real SN = std::sin(kx*x1+ky*x2);
              vxs(tk,tj,ti) = 2.0*vx*vol*SN;
              dvys(tk,tj,ti) = 2.0*dvy*vol*SN;
            }
          }
        }
      } else if (loc.level-1 == root_level) { // level difference 1
        if (pmb->block_size.nx3==1) { // 2D
          int k = pmb->ks;
          for (int j=pmb->cjs; j<=pmb->cje; j++) {
            for (int i=pmb->cis; i<=pmb->cie; i++) {
              int ii = (i-pmb->cis)*2+pmb->is;
              int jj = (j-pmb->cjs)*2+pmb->js;
              int kk = k;
              int ti = static_cast<int>(loc.lx1>>1)*pmb->block_size.nx1
                       +static_cast<int>(loc.lx1%2)*(pmb->block_size.nx1/2)+(i-pmb->cis);
              int tj = static_cast<int>(loc.lx2>>1)*pmb->block_size.nx2
                       +static_cast<int>(loc.lx2%2)*(pmb->block_size.nx2/2)+(j-pmb->cjs);
              int tk = k;
              Real vol00 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii  );
              Real vol01 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii+1);
              Real vol10 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii  );
              Real vol11 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii+1);
              Real vx00 = pmb->phydro->w(IVX,kk  ,jj  ,ii  );
              Real vx01 = pmb->phydro->w(IVX,kk  ,jj  ,ii+1);
              Real vx10 = pmb->phydro->w(IVX,kk  ,jj+1,ii  );
              Real vx11 = pmb->phydro->w(IVX,kk  ,jj+1,ii+1);
              Real dvy00 = pmb->phydro->w(IVY,kk  ,jj  ,ii  );
              Real dvy01 = pmb->phydro->w(IVY,kk  ,jj  ,ii+1);
              Real dvy10 = pmb->phydro->w(IVY,kk  ,jj+1,ii  );
              Real dvy11 = pmb->phydro->w(IVY,kk  ,jj+1,ii+1);
              if(!pmb->porb->orbital_advection_defined) {
                dvy00 += qom*pmb->pcoord->x1v(ii  );
                dvy01 += qom*pmb->pcoord->x1v(ii+1);
                dvy10 += qom*pmb->pcoord->x1v(ii  );
                dvy11 += qom*pmb->pcoord->x1v(ii+1);
              }
              Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                 +ky*pmb->pcoord->x2f(jj+1));
              Real vx_vol  = vx00*vol00+vx01*vol01
                             +vx10*vol10+vx11*vol11;
              Real dvy_vol = dvy00*vol00+dvy01*vol01
                             +dvy10*vol10+dvy11*vol11;
              vxs(tk,tj,ti) = 2.0*SN*vx_vol;
              dvys(tk,tj,ti) = 2.0*SN*dvy_vol;
            }
          }
        } else { // 3D
          for (int k=pmb->cks; k<=pmb->cke; k++) {
            for (int j=pmb->cjs; j<=pmb->cje; j++) {
              for (int i=pmb->cis; i<=pmb->cie; i++) {
                int ii = (i-pmb->cis)*2+pmb->is;
                int jj = (j-pmb->cjs)*2+pmb->js;
                int kk = (k-pmb->cks)*2+pmb->ks;
                int ti = static_cast<int>(loc.lx1>>1)*pmb->block_size.nx1
                         +static_cast<int>(loc.lx1%2)*(pmb->block_size.nx1/2)
                         +(i-pmb->cis);
                int tj = static_cast<int>(loc.lx2>>1)*pmb->block_size.nx2
                         +static_cast<int>(loc.lx2%2)*(pmb->block_size.nx2/2)
                         +(j-pmb->cjs);
                int tk = static_cast<int>(loc.lx3>>1)*pmb->block_size.nx3
                         +static_cast<int>(loc.lx3%2)*(pmb->block_size.nx3/2)
                         +(k-pmb->cks);
                Real vol000 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii  );
                Real vol001 = pmb->pcoord->GetCellVolume(kk  ,jj  ,ii+1);
                Real vol010 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii  );
                Real vol011 = pmb->pcoord->GetCellVolume(kk  ,jj+1,ii+1);
                Real vol100 = pmb->pcoord->GetCellVolume(kk+1,jj  ,ii  );
                Real vol101 = pmb->pcoord->GetCellVolume(kk+1,jj  ,ii+1);
                Real vol110 = pmb->pcoord->GetCellVolume(kk+1,jj+1,ii  );
                Real vol111 = pmb->pcoord->GetCellVolume(kk+1,jj+1,ii+1);
                Real vx000 = pmb->phydro->w(IVX,kk  ,jj  ,ii  );
                Real vx001 = pmb->phydro->w(IVX,kk  ,jj  ,ii+1);
                Real vx010 = pmb->phydro->w(IVX,kk  ,jj+1,ii  );
                Real vx011 = pmb->phydro->w(IVX,kk  ,jj+1,ii+1);
                Real vx100 = pmb->phydro->w(IVX,kk+1,jj  ,ii  );
                Real vx101 = pmb->phydro->w(IVX,kk+1,jj  ,ii+1);
                Real vx110 = pmb->phydro->w(IVX,kk+1,jj+1,ii  );
                Real vx111 = pmb->phydro->w(IVX,kk+1,jj+1,ii+1);
                Real dvy000 = pmb->phydro->w(IVY,kk  ,jj  ,ii  );
                Real dvy001 = pmb->phydro->w(IVY,kk  ,jj  ,ii+1);
                Real dvy010 = pmb->phydro->w(IVY,kk  ,jj+1,ii  );
                Real dvy011 = pmb->phydro->w(IVY,kk  ,jj+1,ii+1);
                Real dvy100 = pmb->phydro->w(IVY,kk+1,jj  ,ii  );
                Real dvy101 = pmb->phydro->w(IVY,kk+1,jj  ,ii+1);
                Real dvy110 = pmb->phydro->w(IVY,kk+1,jj+1,ii  );
                Real dvy111 = pmb->phydro->w(IVY,kk+1,jj+1,ii+1);
                if(!pmb->porb->orbital_advection_defined) {
                  dvy000 += qom*pmb->pcoord->x1v(ii  );
                  dvy001 += qom*pmb->pcoord->x1v(ii+1);
                  dvy010 += qom*pmb->pcoord->x1v(ii  );
                  dvy011 += qom*pmb->pcoord->x1v(ii+1);
                  dvy100 += qom*pmb->pcoord->x1v(ii  );
                  dvy101 += qom*pmb->pcoord->x1v(ii+1);
                  dvy110 += qom*pmb->pcoord->x1v(ii  );
                  dvy111 += qom*pmb->pcoord->x1v(ii+1);
                }
                Real SN = std::sin(kx*pmb->pcoord->x1f(ii+1)
                                   +ky*pmb->pcoord->x2f(jj+1));
                Real vx_vol = vx000*vol000+vx001*vol001
                              +vx010*vol010+vx011*vol011
                              +vx100*vol100+vx101*vol101
                              +vx110*vol110+vx111*vol111;
                Real dvy_vol = dvy000*vol000+dvy001*vol001
                               +dvy010*vol010+dvy011*vol011
                               +dvy100*vol100+dvy101*vol101
                               +dvy110*vol110+dvy111*vol111;
                vxs(tk,tj,ti) = 2.0*SN*vx_vol;
                dvys(tk,tj,ti) = 2.0*SN*dvy_vol;
              }
            }
          }
        }
      } else { // level difference 2
        std::stringstream msg;
        msg << "### FATAL ERROR in ssheet.cpp ProblemGenerator"   << std::endl
            << "This problem prohibits level > 1 for ipert=3 "
            << "with error_output=true."  << std::endl;
          ATHENA_ERROR(msg);
      }
    } // pmb
#ifdef MPI_PARALLEL
    if (Globals::nranks > 1) {
      int ntot = mesh_size.nx3*mesh_size.nx2*mesh_size.nx1;
      // vxs
      if (Globals::my_rank == 0) {
        MPI_Reduce(MPI_IN_PLACE, vxs.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(vxs.data(), vxs.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 0, MPI_COMM_WORLD);
      }
      // dvys
      if (Globals::my_rank == 1) {
        MPI_Reduce(MPI_IN_PLACE, dvys.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 1, MPI_COMM_WORLD);
      } else {
        MPI_Reduce(dvys.data(), dvys.data(), ntot, MPI_ATHENA_REAL,
                   MPI_SUM, 1, MPI_COMM_WORLD);
      }
    }
#endif
  } // flag
  return;
}

void lowerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables for lower side of box
  Real rho0 = pmb->pmy_mesh->ruser_mesh_data[0](6);
  Real s0 = pmb->pmy_mesh->ruser_mesh_data[0](7);
  Real cs = pmb->pmy_mesh->ruser_mesh_data[0](8);
  for (int k = kl; k <= ku; ++k) {
    for (int j = 1; j <= ngh; ++j) {
      for (int i = il; i <= iu; ++i) {
	if (pco->x1v(i) < 0) {
	  // flood lower left of box
	  prim(IDN,k,jl-j,i) = 2 * rho0;
	  prim(IVX,k,jl-j,i) = 0;
	  prim(IVY,k,jl-j,i) = - pco->x1v(i) * s0;
	  prim(IPR,k,jl-j,i) = cs * cs * rho0;
	} else {
	  // do not flood lower right
	  prim(IDN,k,jl-j,i) = prim(IDN,k,jl,i);
	  prim(IVX,k,jl-j,i) = prim(IVX,k,jl,i);
	  prim(IVY,k,jl-j,i) = prim(IVY,k,jl,i);
	  prim(IPR,k,jl-j,i) = prim(IPR,k,jl,i);
	}
      }
    }
  }
}

void upperX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim, FaceField &b, Real time, Real dt, int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  // set primitive variables for upper side of box
  Real rho0 = pmb->pmy_mesh->ruser_mesh_data[0](6);
  Real s0 = pmb->pmy_mesh->ruser_mesh_data[0](7);
  Real cs = pmb->pmy_mesh->ruser_mesh_data[0](8);
  for (int k = kl; k <= ku; ++k) {
    for (int j = 1; j <= ngh; ++j) {
      for (int i = il; i <= iu; ++i) {
        if (pco->x1v(i) >= 0) {
	  // flood upper right of box
	  prim(IDN,k,jl+j,i) = rho0;
	  prim(IVX,k,jl+j,i) = 0;
	  prim(IVY,k,jl+j,i) = - pco->x1v(i) * s0;
	  prim(IPR,k,jl+j,i) = cs * cs * rho0;
        } else {
	  // do not flood upper left
	  prim(IDN,k,jl+j,i) = prim(IDN,k,jl,i);
	  prim(IVX,k,jl+j,i) = prim(IVX,k,jl,i);
	  prim(IVY,k,jl+j,i) = prim(IVY,k,jl,i);
	  prim(IPR,k,jl+j,i) = prim(IPR,k,jl,i);
	}
      }
    }
  }
}


void editArray(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // tester, can user enrolled functions edit the ruser_mesh_data arrays?
  pmb->pmy_mesh->ruser_mesh_data[0](1) = 10;
}


// inactive
void fill(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // custom function to replenish gas at shearing box boundaries
  Real x, y;
  Real l = pmb->pmy_mesh->ruser_mesh_data[0](5);
  Real rho0 = pmb->pmy_mesh->ruser_mesh_data[0](6);
  Real s0 = pmb->pmy_mesh->ruser_mesh_data[0](7);
  Real cs = pmb->pmy_mesh->ruser_mesh_data[0](8);
  Real gamma = pmb->pmy_mesh->ruser_mesh_data[0](9);
  int edgeFlag = -1;

  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      y = pmb->pcoord->x2v(j);
      // is cell close to ymax or ymin?
      if (y < -1 + l) {
        edgeFlag = 0;
      } else if (y > 1 - l) {
        edgeFlag = 1;
      } else {
        continue;
      }
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        x = pmb->pcoord->x1v(i);
        // is cell interior or exterior to perturber?
        if (edgeFlag == 0 && x > 0) {
	  continue; // skip cell, low on right
	} else if (edgeFlag == 1 && x < 0) {
	  continue; // skip cell, high on left
	}
	// cell represents new material, fix
	cons(IM1,k,j,i) = 0;
	cons(IM2,k,j,i) = - s0 * x * 2 * rho0;
	cons(IDN,k,j,i) = 2 * rho0;
	cons(IEN,k,j,i) = cs * cs * 2 * rho0 / (gamma - 1);
      }
    }
  }
}

void sink(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // custom function to remove mass from centre
  Mesh *mesh = pmb->pmy_mesh;
  Real G = mesh->ruser_mesh_data[0](0);
  Real m = mesh->ruser_mesh_data[0](1);
  Real Rs = mesh->ruser_mesh_data[0](2);
  Real gmm = mesh->ruser_mesh_data[0](10);
  Real dx = pmb->pcoord->x1v(1) - pmb->pcoord->x1v(0);
  Real A = dx * dx;
  Real h = mesh->mesh_size.x3max - mesh->mesh_size.x3min; // = 1 by default
  Real V = A * h;
  Real dm = 0;
  Real x1, x2, v1, v2, R, vsqr, vTheta, E;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // do for all particles
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        v1 = prim(IVX,k,j,i);
        v2 = prim(IVY,k,j,i);
        R = sqrt(x1 * x1 + x2 * x2);
        vsqr = v1 * v1 + v2 * v2;
        E = 0.5 * vsqr - G * m / R;
        if (R < Rs && E < 0) {
          // particle within sink AND bound, decrement density and kick
          cons(IDN,k,j,i) /= (1 + gmm * dt);
          vTheta = (x1 * v2 - v1 * x2) / R;
          cons(IM1,k,j,i) += prim(IDN,k,j,i) * gmm * dt * -x2 / R;
          cons(IM2,k,j,i) += prim(IDN,k,j,i) * gmm * dt * x1 / R;
          // accrete to central body
          if (mesh->ruser_mesh_data[0](4) == 0) {
            dm += gmm * dt * cons(IDN,k,j,i) * V;
          }
        }
      }
    }
  }
  // update central body
  pmb->pmy_mesh->ruser_mesh_data[0](1) = 100; // ISSUE IS HERE: cannot seem to edit ruser field
  return;
}

void suck(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // custom function to add central gravity well
  Real x1, x2, R, u, g;
  // grab global data
  Mesh *mesh = pmb->pmy_mesh;
  Real G = mesh->ruser_mesh_data[0](0);
  Real m = mesh->ruser_mesh_data[0](1); // able to load data here, but not pass back?
  Real h = mesh->ruser_mesh_data[0](2);
  Real fac = - G * m / (h * h);
  for (int k = pmb->ks; k<= pmb->ke; ++k) {
    for(int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        // apply central gravitational acceleration to all particles
        x1 = pmb->pcoord->x1v(i);
        x2 = pmb->pcoord->x2v(j);
        R = sqrt(x1 * x1 + x2 * x2);
        u = R / h;
        // determine softening regime
        if (u > 1) {
          // long range, true Newtonian potential
          g = 1 / (u * u);
        } else if (u > 0.5) {
          // intermediate
          g = -0.06667 / (u * u) + 21.33333 * u - 48 * u * u + 38.4 * u * u * u - 10.66667 * u * u * u *u;
        } else {
          // close range
          g = 10.66667 * u - 38.4 * u * u * u + 32 * u * u * u * u;
        }
        g *= fac;
        // apply gravitational acceleration
        cons(IM1,k,j,i) += dt * prim(IDN,k,j,i) * g * x1 / R;
        cons(IM2,k,j,i) += dt * prim(IDN,k,j,i) * g * x2 / R;
      }
    }
  }
  return;
}
/*
void radCool(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim, const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  // custom function to include cooling by radiative diffusion
  Real T, dE, sigma;
  ambT4 = ambT * ambT * ambT * ambT;
  sigma = 5.6704e-8;
  for (int k = pmb->ks; k <= pmb->ke; ++k) {
    for (int j = pmb->js; j <= pmb->je; ++j) {
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        rho = cons(IDN,k,j,i) * rho_unit_;
        egas = cons(IEN,k,j,i) * egas_unit_;
        T = invert(*e_of_rho_T,rho, egas, std::max(es - 1.0, 0.1*es)/3.0, float_1pe*2.0*es/3.0);
        dE = - 2.6667 * sigma * (T * T * T * T) / (kappa * rho); // should be using Sigma, mult by disc thickness?
        cons(IEN,k,j,i) += dE * inv_egas_unit_;
      }
    }
  }
  return;
}
*/
Real CustomTimeStep(MeshBlock *pmb) {
  return pmb->pmy_mesh->ruser_mesh_data[0](3);
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(1);
  return;
}

void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  // assign user output data
  user_out_var(0,ks,js,is) = pmy_mesh->ruser_mesh_data[0](1); // this load here works (with initial mass, but not updated values)
}

namespace {

Real Historydvyc(MeshBlock *pmb, int iout) {
  Real qom = qshear*Omega0;
  Real kx = (TWO_PI/x1size)*(static_cast<Real>(nwx));
  Real ky = (TWO_PI/x2size)*(static_cast<Real>(nwy));
  kx += qom*pmb->pmy_mesh->time*ky;
  Real dvyc = 0.0;
  AthenaArray<Real> volume; // 1D array of volumes
  volume.NewAthenaArray(pmb->ncells1);
  Real tvol = x1size*x2size*x3size;
  for (int k=pmb->ks; k<=pmb->ke  ; k++) {
    for (int j=pmb->js; j<=pmb->je  ; j++) {
      pmb->pcoord->CellVolume(k, j, pmb->is, pmb->ie, volume);
      for (int i=pmb->is; i<=pmb->ie  ; i++) {
        Real x1 = pmb->pcoord->x1v(i);
        Real x2 = pmb->pcoord->x2v(j);
        Real CS = std::cos(kx*x1+ky*x2);
        Real dvy = pmb->phydro->w(IVY,k,j,i);
        if(!pmb->porb->orbital_advection_defined)
          dvy += qom*x1;
        dvyc += volume(i)*2.0*dvy*CS;
      }
    }
  }
  Real dvy0 = iso_cs*amp
              *std::fabs(static_cast<Real>(nwy)/static_cast<Real>(nwx));
  return dvyc/(dvy0*tvol);
}

Real Historyvxs(MeshBlock *pmb, int iout) {
  if (Globals::my_rank != 0) return 0.0;
  if (pmb->lid != 0) return 0.0;
  AthenaArray<Real> &vs = pmb->pmy_mesh->ruser_mesh_data[0];
  Real vxs = 0.0;
  Real vxs_temp;
  int nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int nx2 = pmb->pmy_mesh->mesh_size.nx2;
  int nx3 = pmb->pmy_mesh->mesh_size.nx3;
  Real tvol = x1size*x2size*x3size;
  for (int k=0; k<nx3; k++) {
    for (int i=0; i<nx1; i++) {
      vxs_temp = 0.0;
      for (int j=0; j<nx2; j++) {
        vxs_temp += vs(k,j,i);
      }
      vxs += std::fabs(vxs_temp);
    }
  }
  Real vx0 = iso_cs*amp;
  return vxs/(vx0*tvol);
}

Real Historydvys(MeshBlock *pmb, int iout) {
  int exe_rank_dvy = (Globals::nranks>1)? 1 : 0;
  if (Globals::my_rank != exe_rank_dvy) return 0.0;
  if (pmb->lid != 0) return 0.0;
  AthenaArray<Real> &vs = pmb->pmy_mesh->ruser_mesh_data[1];
  Real dvys = 0.0;
  Real dvys_temp;
  int nx1 = pmb->pmy_mesh->mesh_size.nx1;
  int nx2 = pmb->pmy_mesh->mesh_size.nx2;
  int nx3 = pmb->pmy_mesh->mesh_size.nx3;
  Real tvol = x1size*x2size*x3size;
  for (int k=0; k<nx3; k++) {
    for (int i=0; i<nx1; i++) {
      dvys_temp = 0.0;
      for (int j=0; j<nx2; j++) {
        dvys_temp += vs(k,j,i);
      }
      dvys += std::fabs(dvys_temp);
    }
  }
  Real dvy0 = iso_cs*amp
              *std::fabs(static_cast<Real>(nwy)/static_cast<Real>(nwx));
  return dvys/(dvy0*tvol);
}
} // namespace

//////////////////////////
// ic_curvature.hpp
//////////////////////////
// 
// initial condition generator for gevolution for LTB models
//
// Author: Julian Adamek (Université de Genève & Observatoire de Paris & Queen Mary University of London & Universität Zürich)
//
// Last modified: March 2025
//
//////////////////////////

#ifndef  IC_CURVATURE_HEADER
#define  IC_CURVATURE_HEADER

using namespace LATfield2;

void generateIC_curvature(metadata & sim, icsettings & ic, cosmology & cosmo, const double fourpiG, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_cdm, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_b, Particles<part_simple,part_simple_info,part_simple_dataType> * pcls_ncdm, double * maxvel, Field<Real> * phi, Field<Real> * chi, Field<Real> * Bi, Field<Real> * source, Field<Real> * Sij, Field<Cplx> * zetaFT, Field<Cplx> * scalarFT, Field<Cplx> * BiFT, Field<Cplx> * SijFT, PlanFFT<Cplx> * plan_phi, PlanFFT<Cplx> * plan_chi, PlanFFT<Cplx> * plan_Bi, PlanFFT<Cplx> * plan_source, PlanFFT<Cplx> * plan_Sij, parameter * params, int & numparam)
{
	double a = 1. / (1. + sim.z_in);
	float * pcldata = NULL;
	double * sinc = NULL;
	Site x(phi->lattice());
	rKSite kFT(scalarFT->lattice());
	double max_displacement;
	double r2, d1, i2, inner_radius, H2, dlnm;
	part_simple_info pcls_cdm_info;
	part_simple_dataType pcls_cdm_dataType;
	part_simple_info pcls_b_info;
	part_simple_dataType pcls_b_dataType;
	Real boxSize[3] = {1.,1.,1.};
	Real vertex[3];
	Real origin[3] = {Real(0.5 * sim.numpts), Real(0.5 * sim.numpts), Real(0.5 * sim.numpts)};
	Field<Real> * ic_fields[2];
	char filename[2*PARAM_MAX_LENGTH+64];
	int reduce = MAX;
	
	ic_fields[0] = chi;
	ic_fields[1] = chi;

	H2 = Hconf(a, fourpiG, cosmo) * Hconf(a, fourpiG, cosmo);
	d1 = -0.6 * ic.LTB_Omega_k * (1. + 11. * ic.LTB_Omega_k / 35.);
	inner_radius = sim.LTB_radius / pow(1. + d1, 1./3.);
	i2 = inner_radius * inner_radius;
	dlnm = 0.5 * d1 * H2 * sim.LTB_radius * sim.LTB_radius;

	COUT << " computed gravitational potential at the center of the LTB model: phi(r=0) = " << -0.25 * d1 * H2 * (sim.LTB_radius * sim.LTB_radius * (1.0 - d1 / 3.)) << endl; //<< -0.5 * dlnm << endl;

	dlnm *= 4. * M_PI * sim.LTB_radius * sim.LTB_radius * sim.LTB_radius / 3.;
	
	// generate the kernel for the 1st order displacement field

	loadHomogeneousTemplate(ic.pclfile[0], sim.numpcl[0], pcldata);
	
	if (pcldata == NULL)
	{
		COUT << " error: particle data was empty!" << endl;
		parallel.abortForce();
	}
	
	if (ic.flags & ICFLAG_CORRECT_DISPLACEMENT)
		generateCICKernel(*source, sim.numpcl[0], pcldata, ic.numtile[0]);
	else
		generateCICKernel(*source);	
	
	plan_source->execute(FFT_FORWARD);
	

	// use BiFT as temporary storage for the kernel
	for (kFT.first(); kFT.test(); kFT.next())
		(*BiFT)(kFT, 0) = (*scalarFT)(kFT);
	
	// precompute sinc function
	sinc = (double *) malloc(sim.numpts * sizeof(double));

	sinc[0] = 1.;
	for (int i = 1; i <= sim.numpts/2; i++)
		sinc[i] = sin(M_PI * i / sim.numpts) / (M_PI * i / sim.numpts);

	for (int i = 1; i < sim.numpts/2; i++)
		sinc[sim.numpts-i] = sinc[i];

	for (x.first(); x.test(); x.next())
	{
		r2 = 0.;
		for (int i = 0; i < 3; i++)
		{
			r2 += (x.coord(i) - sim.numpts / 2) * (x.coord(i) - sim.numpts / 2);
			vertex[i] = (Real) x.coord(i);
		}

		r2 /= sim.numpts * sim.numpts;

		Real w1 = Real(1) - computeTruncatedCellVolume(vertex, origin, inner_radius * sim.numpts);
		Real w2 = Real(1) - computeTruncatedCellVolume(vertex, origin, sim.LTB_radius * sim.numpts);

		/*if (w1 > 1 || w1 < -1.0e-5)
			cout << "error in position [" << x.coord(0) << "," << x.coord(1) << "," << x.coord(2) << "]: w1 = " << w1 << endl;
		if (w2 > 1 || w2 < -1.0e-5)
			cout << "error in position [" << x.coord(0) << "," << x.coord(1) << "," << x.coord(2) << "]: w2 = " << w2 << endl;*/

		w2 -= w1;
		
		(*source)(x) = w1 * d1 * (1. + 0.5 * sim.gr_flag * H2 * (sim.LTB_radius * sim.LTB_radius * (1. + 10. * d1 / 21.) - r2 * (1. + 16. * d1 / 63.) + d1 * H2 * (17. * r2 * r2 + sim.LTB_radius * sim.LTB_radius * (27. * sim.LTB_radius * sim.LTB_radius - 44. * r2)) / 24.)) - w2;
		
		if (r2 < i2)
		{
			(*phi)(x) = -0.25 * d1 * H2 * (sim.LTB_radius * sim.LTB_radius * (1.0 - d1 / 3.) - r2 - 5. * d1 * H2 * r2 * (2. * r2 + sim.LTB_radius * sim.LTB_radius - 3. * sim.LTB_radius * sqrt(r2)) / 6.);
		}
		else if (r2 < sim.LTB_radius * sim.LTB_radius)
		{
			r2 = sqrt(r2);
			(*phi)(x) = -0.25 * H2 * (1. + 2. * sim.LTB_radius / r2) * (sim.LTB_radius - r2) * (sim.LTB_radius - r2);
		}
		else
		{
			(*phi)(x) = Real(0);
		}
	}

	plan_phi->execute(FFT_FORWARD);

	for (kFT.first(); kFT.test(); kFT.next())
	{
		double W = sinc[kFT.coord(0)] * sinc[kFT.coord(1)] * sinc[kFT.coord(2)];
		(*scalarFT)(kFT) *= W*W/sim.numpts/sim.numpts/sim.numpts;
	}

	plan_phi->execute(FFT_BACKWARD);
	phi->updateHalo();
	
	sprintf(filename, "%s%s_IC_phi.h5", sim.output_path, sim.basename_generic);
	phi->saveHDF5(string(filename));
	
    source->updateHalo();

	sprintf(filename, "%s%s_IC_delta.h5", sim.output_path, sim.basename_generic);
	source->saveHDF5(string(filename));

	COUT << " computing 1st order displacement..." << endl;
	
	if (sim.gr_flag == 1)
	{
		for (x.first(); x.test(); x.next())
			(*chi)(x) = (*source)(x) - 3. * (*phi)(x);
	}
	else
	{
		for (x.first(); x.test(); x.next())
			(*chi)(x) = (*source)(x);
	}

	plan_chi->execute(FFT_FORWARD);
	
	for (kFT.first(); kFT.test(); kFT.next())
	{
		if ((*BiFT)(kFT, 0).norm() > 1.0e-16)
			(*scalarFT)(kFT) = (*scalarFT)(kFT) / (*BiFT)(kFT, 0);
	}
	
	plan_chi->execute(FFT_BACKWARD);
	chi->updateHalo();
	
	//filename.assign("1st_order_displacement.h5");
	sprintf(filename, "%s%s_IC_1st_order_displacement.h5", sim.output_path, sim.basename_generic);
	chi->saveHDF5(string(filename));
	
	strcpy(pcls_cdm_info.type_name, "part_simple");
	if (sim.baryon_flag == 1)
		pcls_cdm_info.mass = (1. + dlnm) * cosmo.Omega_cdm / (Real) (sim.numpcl[0]*(long)ic.numtile[0]*(long)ic.numtile[0]*(long)ic.numtile[0]);
	else
		pcls_cdm_info.mass = (1. + dlnm) * (cosmo.Omega_cdm + cosmo.Omega_b) / (Real) (sim.numpcl[0]*(long)ic.numtile[0]*(long)ic.numtile[0]*(long)ic.numtile[0]);
	pcls_cdm_info.relativistic = false;
	
	pcls_cdm->initialize(pcls_cdm_info, pcls_cdm_dataType, &(phi->lattice()), boxSize);
	
	initializeParticlePositions(sim.numpcl[0], pcldata, ic.numtile[0], *pcls_cdm);
	
	free(pcldata);

	if (sim.baryon_flag == 3)	// baryon treatment = hybrid; displace particles using both displacement fields
		pcls_cdm->moveParticles(displace_pcls_ic_basic, 1./sim.numpts/sim.numpts/sim.numpts, ic_fields, 2, NULL, &max_displacement, &reduce, 1);
	else
		pcls_cdm->moveParticles(displace_pcls_ic_basic, 1./sim.numpts/sim.numpts/sim.numpts, &chi, 1, NULL, &max_displacement, &reduce, 1);	// displace CDM particles
	
	sim.numpcl[0] *= (long) ic.numtile[0] * (long) ic.numtile[0] * (long) ic.numtile[0];
	
	COUT << " " << sim.numpcl[0] << " cdm particles initialized: maximum displacement (1st order) = " << max_displacement * sim.numpts << " lattice units." << endl;

	maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, a/(1.5 * Hconf(a, fourpiG, cosmo)), &phi, 1) / a;

	parallel.max<double>(maxvel, 1);

	COUT << " maximum velocity (1st order) = " << maxvel[0] << endl;

	COUT << " computing 1st order density..." << endl;
	
	projection_init(chi);
	projection_T00_project(pcls_cdm, chi, a, sim.gr_flag ? phi : NULL); // gr_flag=1 -> GR, else Newton
	scalarProjectionCIC_comm(chi);
	
	for (x.first(); x.test(); x.next())
		(*chi)(x) = (*chi)(x) / (cosmo.Omega_cdm + cosmo.Omega_b) - 1.;
		
	//filename.assign("1st_order_density.h5");
	sprintf(filename, "%s%s_IC_1st_order_density.h5", sim.output_path, sim.basename_generic);
	chi->saveHDF5(string(filename));

	for (x.first(); x.test(); x.next())
		(*chi)(x) = (*source)(x) - (*chi)(x);
		
	plan_chi->execute(FFT_FORWARD);
	
	for (kFT.first(); kFT.test(); kFT.next())
	{
		if ((*BiFT)(kFT, 0).norm() > 1.0e-16)
			(*scalarFT)(kFT) = (*scalarFT)(kFT) / (*BiFT)(kFT, 0);
	}
	
	plan_chi->execute(FFT_BACKWARD);
	chi->updateHalo();
	
	//filename.assign("2nd_order_displacement.h5");
	//chi->saveHDF5(filename);
	
	pcls_cdm->moveParticles(displace_pcls_ic_basic, 1./sim.numpts/sim.numpts/sim.numpts, &chi, 1, NULL, &max_displacement, &reduce, 1);
	
	COUT << " " << sim.numpcl[0] << " cdm particles initialized: maximum displacement (2nd order) = " << max_displacement * sim.numpts << " lattice units." << endl;

	// compute velocity potential at second order
	for (x.first(); x.test(); x.next())
	{
		r2 = 0.;
		for (int i = 0; i < 3; i++)
			r2 += (x.coord(i) - sim.numpts / 2) * (x.coord(i) - sim.numpts / 2);

		r2 /= sim.numpts * sim.numpts;

		(*chi)(x) = (*phi)(x) * (Real(1) + Real(2) * (*phi)(x));

		if (r2 < sim.LTB_radius * sim.LTB_radius)
			(*chi)(x) -= d1 * d1 * H2 * (r2 - sim.LTB_radius * sim.LTB_radius) * (16 + 35 * H2 * (r2 - sim.LTB_radius * sim.LTB_radius)) / 336.;
	}

	chi->updateHalo();

	maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, a/(1.5 * Hconf(a, fourpiG, cosmo)), &chi, 1) / a;

	parallel.max<double>(maxvel, 1);

	COUT << " maximum velocity (2nd order) = " << maxvel[0] << endl;

	COUT << " computing 2nd order density..." << endl;

    projection_init(chi);
    projection_T00_project(pcls_cdm, chi, a, sim.gr_flag ? phi : NULL); // toma
    scalarProjectionCIC_comm(chi);

    for (x.first(); x.test(); x.next())
            (*chi)(x) = (*source)(x) - ((*chi)(x) / (cosmo.Omega_cdm + cosmo.Omega_b) - 1.);

	plan_chi->execute(FFT_FORWARD);

    for (kFT.first(); kFT.test(); kFT.next())
    {
        if ((*BiFT)(kFT, 0).norm() > 1.0e-16)
            (*scalarFT)(kFT) = (*scalarFT)(kFT) / (*BiFT)(kFT, 0);
    }

    plan_chi->execute(FFT_BACKWARD);
    chi->updateHalo();

    //filename.assign("3rd_order_displacement.h5");
    //chi->saveHDF5(filename);

	pcls_cdm->moveParticles(displace_pcls_ic_basic, 1./sim.numpts/sim.numpts/sim.numpts, &chi, 1, NULL, &max_displacement, &reduce, 1);

    COUT << " " << sim.numpcl[0] << " cdm particles initialized: maximum displacement (3rd order) = " << max_displacement * sim.numpts << " lattice units." << endl;
	
	/*if (sim.baryon_flag == 1)
	{
		loadHomogeneousTemplate(ic.pclfile[1], sim.numpcl[1], pcldata);
	
		if (pcldata == NULL)
		{
			COUT << " error: particle data was empty!" << endl;
			parallel.abortForce();
		}
		
		strcpy(pcls_b_info.type_name, "part_simple");
		pcls_b_info.mass = cosmo.Omega_b / (Real) (sim.numpcl[1]*(long)ic.numtile[1]*(long)ic.numtile[1]*(long)ic.numtile[1]);
		pcls_b_info.relativistic = false;
	
		pcls_b->initialize(pcls_b_info, pcls_b_dataType, &(phi->lattice()), boxSize);
	
		initializeParticlePositions(sim.numpcl[1], pcldata, ic.numtile[1], *pcls_b);
		
		pcls_b->moveParticles(displace_pcls_ic_basic, 1./sim.boxsize/sim.boxsize, &phi, 1, NULL, &max_displacement, &reduce, 1);	// displace baryon particles
	
		sim.numpcl[1] *= (long) ic.numtile[1] * (long) ic.numtile[1] * (long) ic.numtile[1];
	
		COUT << " " << sim.numpcl[1] << " baryon particles initialized: maximum displacement = " << max_displacement * sim.numpts << " lattice units." << endl;
	
		free(pcldata);
	}
	
	if (ic.pkfile[0] == '\0')	// set velocities using transfer functions
	{
		filename.assign(ic.velocityfile[1]);
		
		chi->loadHDF5(filename);
		
		chi->updateHalo();
		
		if (sim.baryon_flag == 3)	// baryon treatment = hybrid; set velocities using both velocity potentials
			maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, -a/sim.boxsize, ic_fields, 2) / a;
		else
			maxvel[0] = pcls_cdm->updateVel(initialize_q_ic_basic, -a/sim.boxsize, &chi, 1) / a;	// set CDM velocities
		
		if (sim.baryon_flag == 1)
			maxvel[1] = pcls_b->updateVel(initialize_q_ic_basic, -a/sim.boxsize, &phi, 1) / a;	// set baryon velocities
	}*/
	
	if (sim.baryon_flag > 1) sim.baryon_flag = 0;
	
	projection_init(Bi);
	projection_T0i_project(pcls_cdm, Bi,  phi);
	if (sim.baryon_flag)
		projection_T0i_project(pcls_b, Bi, phi);
	projection_T0i_comm(Bi);
	prepareFTsource(*Bi, *phi, 3. * a * a * Hconf(a, fourpiG, cosmo) * (double) sim.numpts / fourpiG);
	plan_Bi->execute(FFT_FORWARD);
	projectFTvector(*BiFT, *BiFT, fourpiG / (double) sim.numpts / (double) sim.numpts);	
	plan_Bi->execute(FFT_BACKWARD);	
	Bi->updateHalo();	// B initialized
	
	ic_fields[1] = Bi;
	
	projection_init(Sij);
	projection_Tij_project(pcls_cdm, Sij, a, phi);
	if (sim.baryon_flag)
		projection_Tij_project(pcls_b, Sij, a, phi);
	projection_Tij_comm(Sij);
	
	prepareFTsource<Real>(*phi, *Sij, *Sij, 2. * fourpiG / a / (double) sim.numpts / (double) sim.numpts);	
	plan_Sij->execute(FFT_FORWARD);	
	projectFTscalar(*SijFT, *scalarFT);
	plan_chi->execute(FFT_BACKWARD);		
	chi->updateHalo();	// chi now finally contains chi

	free(sinc);
}

#endif

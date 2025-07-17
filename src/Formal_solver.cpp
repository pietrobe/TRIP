#include "Formal_solver.hpp"

/////////////////////// for BESSER

inline Real PSI_M_LIN(Real ex, Real t) {
    return ( t > 0.11 ? ((1.0-ex*(1.0+t))/t) : ((t*(t*(t*(t*(t*(t*((63.0-8.0*t)*t-432.0)+2520.0)-12096.0)+45360.0)-120960.0)+181440.0))/362880.0) );
}

inline Real PSI_O_LIN(Real ex, Real t) {
    return ( t > 0.11 ? ((ex+t-1.0)/t)       : ((t*(t*(t*(t*(t*(t*((9.0-t)*t-72.0)+504.0)-3024.0)+15120.0)-60480.0)+181440.0))/362880.0) );
}

inline Real OMEGA_M(Real ex, Real t) {
    return ( t > 0.14 ? ((2.0-ex*(t*t+2.0*t+2.0))/(t*t))    : ((t*(t*(t*(t*(t*(t*((140.0-18.0*t)*t-945.0)+5400.0)-25200.0)+90720.0)-226800.0)+302400.0))/907200.0) );
}

inline Real OMEGA_O(Real ex, Real t) {
    return ( t > 0.18 ? (1.0-2.0*(t+ex-1.0)/(t*t))          : ((t*(t*(t*(t*(t*(t*((10.0-t)*t-90.0)+720.0)-5040.0)+30240.0)-151200.0)+604800.0))/1814400.0) );
}

inline Real OMEGA_C(Real ex, Real t) {
    return ( t > 0.18 ? (2.0*(t-2.0+ex*(t+2.0))/(t*t))      : ((t*(t*(t*(t*(t*(t*((35.0-4.0*t)*t-270.0)+1800.0)-10080.0)+45360.0)-151200.0)+302400.0))/907200.0) );
}

inline int YBETWAB(Real y, Real a, Real b) {
    return (((a<=b && y>=a && y<=b) || (a>=b && y<=a && y>=b)) ? 1 : 0);
}

// some basic scalar functions for the SC method (Kunasz Auer 1988)
inline Real SC_linear(const Real dt, const Real I_in, const Real S1, const Real S2)
{	
	const long double dt_aux = -dt;     
	const long double e_dt = std::exp(-dt_aux); // TOD check this
	const long double u0 = 1.0 - e_dt;
	const long double u1 = dt_aux - u0;

	const long double u1_dt_m = u1 / dt_aux;
	    
	return e_dt * I_in + S1 * (u0 - u1_dt_m) + S2 * u1_dt_m;  
}


inline Real SC_parabolic(const Real dt1, const Real dt2, const Real I_in, const Real S1, const Real S2, const Real S3)
{
	const long double dt_aux1 = -dt1;     
	const long double dt_aux2 = -dt2;

	const long double e_dt = std::exp(-dt_aux1);
	const long double u0 = 1.0 - e_dt;
	const long double u1 = dt_aux1 - u0;
    const long double u2 = dt_aux1 * dt_aux1 - 2 * u1;              
       
    const long double phi1 = u0 + (u2 - (dt_aux2 + 2 * dt_aux1) * u1)/(dt_aux1 * (dt_aux1 + dt_aux2));
    const long double phi2 = ((dt_aux1 + dt_aux2) * u1 - u2) / (dt_aux1 * dt_aux2);
    const long double phi3 = (u2 - dt_aux1 * u1) / (dt_aux2 * (dt_aux1 + dt_aux2));
   
    return e_dt * I_in + S1 * phi1 + S2 * phi2 + S3 * phi3;  
}


Real CorrectYAB(Real y, Real a, Real b) {
    Real min = fmin(a,b), max = fmax(a,b);
    if (y < min) return min;
    else if (y > max) return max;
    else return y;
}

Real mat_QBezierC0(Real h0, Real h1, Real ym, Real yo, Real yp) {
    Real dm, dp, yder, c0, c1;
    int cond0, cond1;

    if (h0>0.0 && h1>0.0) {
        dm = (yo - ym) / h0;
        dp = (yp - yo) / h1;
    }
    else {
        c0 = yo;
        return c0;
    }

    if (dm * dp <= 0) {
        c0 = yo;
        return c0;
    }

    yder = (h0*dp + h1*dm)/(h0 + h1);
    c0 = yo - 0.5 * h0 * yder;
    c1 = yo + 0.5 * h1 * yder;

    cond0 = YBETWAB(c0,ym,yo);
    cond1 = YBETWAB(c1,yo,yp);
    if (cond0 && cond1) {
        return c0;
    }
    else if (!cond0) {
        c0 = CorrectYAB(c0, ym, yo);
    }
    else if (!cond1) {
        c1 = CorrectYAB(c1, yo, yp);
        yder = 2.0*(c1-yo)/h1;
        c0 = yo - 0.5 * h0 * yder;
        if (!YBETWAB(c0, ym, yo)) {
            c0 = CorrectYAB(c0, ym, yo);
        }
    }
    else {
        std::cerr << "ERROR mat_QBezierC0!" << std::endl;
    }

    if (!YBETWAB(c0, ym, yo)) {
        std::cerr << "ERROR mat_QBezierC0!" << std::endl;
    }

    return c0;
}
///////////////////////

// solve I' = KI - S

void Formal_solver::one_step(const Real dt, input_vec &K1, input_vec &K2, input_vec &S1, input_vec &S2, input_vec &I_in, output_vec &I_out){

	// sanity checks
	if (debug_mode_ and mpi_rank_ == 0)
	{
		if (K2.size() != 16 or S2.size() != 4 or I_in.size() != 4)  std::cerr << "ERROR: wrong input size in one_step().\n";		
	}
	
	std::vector<Real> A(16); 
	std::vector<Real> b(4); 
	
	if (type_ == "implicit_Euler") // K1 and S1 are not needed here 
	{		
		// compute b = I_in - dt * S and A = Id - dt * K		
		int index_ij;

		for (int i = 0; i < 4; ++i) 
		{
			b[i] = I_in[i] - dt * S2[i];						

			for (int j = 0; j < 4; ++j)
			{
				index_ij = 4 * i + j;

				A[index_ij] = (i == j) - dt * K2[index_ij];
			}			
		}		
			
		I_out = solve_4_by_4_system(A, b);			
	}
	else if (type_ == "trapezoidal" or type_ == "Crank–Nicolson")
	{		
		if (debug_mode_ and (K1.size() != 16 or S1.size() != 4) and mpi_rank_ == 0) std::cerr << "\nERROR: wrong input size for Crank–Nicolson!\n";

		std::vector<Real> K_times_I_in(4); 

		Real dt_half = 0.5 * dt;

		int index_ij;

		for (int i = 0; i < 4; ++i) 			
		{						
			for (int j = 0; j < 4; ++j)
			{
				index_ij = 4 * i + j;

				// Id - (dt/2) * K
				A[index_ij] = (i == j) - dt_half * K2[index_ij];

				// K * I_in
				K_times_I_in[i] += K1[index_ij] * I_in[j]; 
			}	

			b[i] = I_in[i] + dt_half * (K_times_I_in[i] - S1[i] - S2[i]);			
		}
			
		I_out = solve_4_by_4_system(A, b);	
	}
	else // DELO_linear is used if not specified
	{    			
		if (debug_mode_ and (K1.size() != 16 or S1.size() != 4) and mpi_rank_ == 0) std::cerr << "\nERROR: wrong input size for DELO_linear!\n";
		
		// solve I' = - KI + S
		const long double dt_aux = -dt;         
        const long double E = std::exp(-dt_aux);      
        const long double F = 1.0 - E;
        const long double G = (1.0 - (1.0 + dt_aux) * E) / dt_aux; 
        
        std::vector<Real> K_times_I_in(4, 0.0);
        int index_ij;
				
        bool Id;

     	for (int i = 0; i < 4; ++i) 
		{			
			for (int j = 0; j < 4; ++j)
			{
				index_ij = 4 * i + j;

				Id = (i == j);
			
				A[index_ij] = Id + (F - G) * (K2[index_ij] - Id);

				// K * I_in
				K_times_I_in[i] += (K1[index_ij] - Id) * I_in[j]; 								
			}			

			b[i] = E * I_in[i] - G * K_times_I_in[i] + G * S1[i] + (F - G) * S2[i];
		}

		I_out = solve_4_by_4_system(A, b);	
		// I_out = solve_4_by_4_system_optimized(A,b);
	}	
	
	// // check valore negativo
	// if (debug_mode_ and I_out[0] < 0 and I_in[0] >= 0 and S1[0] > 0  and S2[0] > 0) std::cout << "\nWARNING: negative intensity detected!\n";	

	// if (I_out[0] < 0) I_out[0] = 0;		
}

// scalar equations
Real Formal_solver::one_step(const Real dt, const Real I_in, const Real S1, const Real S2)
{
	// for now only SC_linear 	
	return SC_linear(dt, I_in, S1, S2);
}


Real Formal_solver::one_step_quadratic(const Real dt1, const Real dt2, const Real I_in, const Real S1, const Real S2, const Real S3)
{
	// for now only SC_quadratic 
	return SC_parabolic(dt1, dt2, I_in, S1, S2, S3);
}


void Formal_solver::one_step_quadratic(const Real dt_1, const Real dt_2, input_vec &K1, input_vec &K2, input_vec &K3,
								  							   input_vec &S1, input_vec &S2, input_vec &S3,
								 							   input_vec &I_in, output_vec &I_out){
	// sanity checks
	if (debug_mode_ and mpi_rank_ == 0)
	{
		if (K1.size() != 16 or K2.size() != 16 or K3.size() != 16) std::cerr << "ERROR: wrong input size of K in one_step().\n"; 
		if (S1.size() != 4  or S2.size() != 4  or S3.size() != 4)  std::cerr << "ERROR: wrong input size of S in one_step().\n";
		if (I_in.size() != 4)   								   std::cerr << "ERROR: wrong input size of I in one_step().\n";
		
		if (stencil_size_ != 3) std::cerr << "ERROR: stencil_size_ should be 3 in one_step().\n";
	}

	if (type_ == "BESSER")
	{
		static const Real VACUUM_OPACITY = 1e-30;
	    static const int ETA_I = 0, ETA_Q = 1, ETA_U = 2, ETA_V = 3, RHO_Q = 11, RHO_U = 13, RHO_V = 6;
	    static const int STOKES_I = 0, STOKES_Q = 1, STOKES_U = 2, STOKES_V = 3;
	 
	    Real tm = - dt_1 + VACUUM_OPACITY, tp = - dt_2 + VACUUM_OPACITY; // added minus since input dt are negative
	    Real ex = exp(-tm);
	    Real om_m = OMEGA_M(ex, tm), om_o = OMEGA_O(ex, tm), om_c = OMEGA_C(ex, tm);
	    Real vec[4], kappa[4][4], id, c0[4];
	    int i, j;
	    Real psi_m_lin = PSI_M_LIN(ex,tm), psi_o_lin = PSI_O_LIN(ex,tm);
	    Real a = -psi_m_lin * K1[ETA_Q]/(K1[ETA_I]+VACUUM_OPACITY);
	    Real b = -psi_m_lin * K1[ETA_U]/(K1[ETA_I]+VACUUM_OPACITY);
	    Real c = -psi_m_lin * K1[ETA_V]/(K1[ETA_I]+VACUUM_OPACITY);
	    Real s = -psi_m_lin * K1[RHO_V]/(K1[ETA_I]+VACUUM_OPACITY);
	    Real q = -psi_m_lin * K1[RHO_U]/(K1[ETA_I]+VACUUM_OPACITY);
	    Real r = -psi_m_lin * K1[RHO_Q]/(K1[ETA_I]+VACUUM_OPACITY);

	    c0[0] = mat_QBezierC0(tm, tp, S1[STOKES_I], S2[STOKES_I], S3[STOKES_I]);
	    c0[1] = mat_QBezierC0(tm, tp, S1[STOKES_Q], S2[STOKES_Q], S3[STOKES_Q]);
	    c0[2] = mat_QBezierC0(tm, tp, S1[STOKES_U], S2[STOKES_U], S3[STOKES_U]);
	    c0[3] = mat_QBezierC0(tm, tp, S1[STOKES_V], S2[STOKES_V], S3[STOKES_V]);

	    // // test to recover DELO linear
	    // c0[0] = 0.5 * (S1[STOKES_I] + S2[STOKES_I]);
		// c0[1] = 0.5 * (S1[STOKES_Q] + S2[STOKES_Q]);
	    // c0[2] = 0.5 * (S1[STOKES_U] + S2[STOKES_U]);
	    // c0[3] = 0.5 * (S1[STOKES_V] + S2[STOKES_V]);

	    vec[0] = ex * I_in[STOKES_I]    +                      a * I_in[STOKES_Q] + b * I_in[STOKES_U] + c * I_in[STOKES_V];
	    vec[1] = ex * I_in[STOKES_Q]    + a * I_in[STOKES_I]                      + s * I_in[STOKES_U] - q * I_in[STOKES_V];
	    vec[2] = ex * I_in[STOKES_U]    + b * I_in[STOKES_I] - s * I_in[STOKES_Q]                      + r * I_in[STOKES_V];
	    vec[3] = ex * I_in[STOKES_V]    + c * I_in[STOKES_I] + q * I_in[STOKES_Q] - r * I_in[STOKES_U];

	    vec[0] += om_m * S1[STOKES_I] + om_o * S2[STOKES_I] + om_c * c0[0];
	    vec[1] += om_m * S1[STOKES_Q] + om_o * S2[STOKES_Q] + om_c * c0[1];
	    vec[2] += om_m * S1[STOKES_U] + om_o * S2[STOKES_U] + om_c * c0[2];
	    vec[3] += om_m * S1[STOKES_V] + om_o * S2[STOKES_V] + om_c * c0[3];

	    id = psi_o_lin /  (K2[ETA_I]+VACUUM_OPACITY);
	    a = id * K2[ETA_Q];
	    b = id * K2[ETA_U];
	    c = id * K2[ETA_V];
	    s = id * K2[RHO_V];
	    q = id * K2[RHO_U];
	    r = id * K2[RHO_Q];

	    id = -c*c*s*s+s*s-2.0*a*c*r*s-2.0*b*c*q*s-a*a*r*r+r*r-2.0*a*b*q*r-b*b*q*q+q*q-c*c-b*b-a*a+1.0;

	    if (0.0 == id) std::cerr << "ERROR in one_step!" << std::endl;
	   
	    id = 1.0/id;

	    kappa[0][0] = s*s+r*r+q*q+1.0;                  kappa[0][1] = -c*r*s-b*s-a*r*r-b*q*r+c*q-a;   kappa[0][2] = -c*q*s+a*s-a*q*r-c*r-b*q*q-b;      kappa[0][3] = -c*s*s-a*r*s-b*q*s+b*r-a*q-c;
	    kappa[1][0] = -c*r*s+b*s-a*r*r-b*q*r-c*q-a;     kappa[1][1] = r*r-c*c-b*b+1.0;                kappa[1][2] = c*c*s-s+q*r+a*c*r+b*c*q+a*b;       kappa[1][3] = r*s-b*c*s-a*b*r-b*b*q+q+a*c;
	    kappa[2][0] = -c*q*s-a*s-a*q*r+c*r-b*q*q-b;     kappa[2][1] = -c*c*s+s+q*r-a*c*r-b*c*q+a*b;   kappa[2][2] = q*q-c*c-a*a+1.0;                   kappa[2][3] = q*s+a*c*s+a*a*r-r+a*b*q+b*c;
	    kappa[3][0] = -c*s*s-a*r*s-b*q*s-b*r+a*q-c;     kappa[3][1] = r*s+b*c*s+a*b*r+b*b*q-q+a*c;    kappa[3][2] = q*s-a*c*s-a*a*r+r-a*b*q+b*c;       kappa[3][3] = s*s-b*b-a*a+1.0;
	    
	    for (j=0; j<4; j++) for (i=0; i<4; i++) kappa[j][i] *= id;

	    I_out[STOKES_I] = kappa[0][0]*vec[0] + kappa[0][1]*vec[1] + kappa[0][2]*vec[2] + kappa[0][3]*vec[3];
	    I_out[STOKES_Q] = kappa[1][0]*vec[0] + kappa[1][1]*vec[1] + kappa[1][2]*vec[2] + kappa[1][3]*vec[3];
	    I_out[STOKES_U] = kappa[2][0]*vec[0] + kappa[2][1]*vec[1] + kappa[2][2]*vec[2] + kappa[2][3]*vec[3];
	    I_out[STOKES_V] = kappa[3][0]*vec[0] + kappa[3][1]*vec[1] + kappa[3][2]*vec[2] + kappa[3][3]*vec[3];
	}
	else if (debug_mode_ and mpi_rank_ == 0)
	{
		std::cerr << "ERROR: formal_solver not supported in one_step_quadratic()\n";
	}
}


void Formal_solver::solve(input_vec &dts, input_field &K, input_field &S, input_vec &I_in, output_field &I_out){

	size_t N_tau = dts.size() + 1;

	// check input sizes
	if (debug_mode_ and (K.size() != N_tau or S.size() != N_tau or I_in.size() != 4) and mpi_rank_ == 0) std::cerr << "\nERROR: wrong input size in Formal_solver::solve()!\n";

	I_out.resize(N_tau);

	// set initial condition 
	I_out[0] = I_in;

	// solve ODE
	for (size_t i = 1; i < N_tau; ++i)
	{
		one_step(dts[i-1], K[i-1], K[i], S[i-1], S[i], I_out[i-1], I_out[i]);		
	}
}












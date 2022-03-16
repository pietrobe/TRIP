#include "Utilities.hpp"

int Mod(int a, int b) 
{
  int r = a % b;
  return r < 0 ? r + b : r;
}

//////////////////////////////// used in W3JS 

static double g_fact[171] = {
  1.000000000000000e+00,  1.000000000000000e+00,  2.000000000000000e+00,  6.000000000000000e+00,
  2.400000000000000e+01,  1.200000000000000e+02,  7.200000000000000e+02,  5.040000000000000e+03,
  4.032000000000000e+04,  3.628800000000000e+05,  3.628800000000000e+06,  3.991680000000000e+07,
  4.790016000000000e+08,  6.227020800000000e+09,  8.717829120000000e+10,  1.307674368000000e+12,
  2.092278988800000e+13,  3.556874280960000e+14,  6.402373705728000e+15,  1.216451004088320e+17,
  2.432902008176640e+18,  5.109094217170944e+19,  1.124000727777608e+21,  2.585201673888498e+22,
  6.204484017332394e+23,  1.551121004333098e+25,  4.032914611266057e+26,  1.088886945041835e+28,
  3.048883446117138e+29,  8.841761993739701e+30,  2.652528598121911e+32,  8.222838654177924e+33,
  2.631308369336936e+35,  8.683317618811887e+36,  2.952327990396042e+38,  1.033314796638615e+40,
  3.719933267899012e+41,  1.376375309122635e+43,  5.230226174666011e+44,  2.039788208119744e+46,
  8.159152832478980e+47,  3.345252661316382e+49,  1.405006117752880e+51,  6.041526306337383e+52,
  2.658271574788449e+54,  1.196222208654802e+56,  5.502622159812090e+57,  2.586232415111682e+59,
  1.241391559253607e+61,  6.082818640342675e+62,  3.041409320171338e+64,  1.551118753287382e+66,
  8.065817517094388e+67,  4.274883284060025e+69,  2.308436973392414e+71,  1.269640335365828e+73,
  7.109985878048636e+74,  4.052691950487721e+76,  2.350561331282878e+78,  1.386831185456898e+80,
  8.320987112741390e+81,  5.075802138772248e+83,  3.146997326038794e+85,  1.982608315404440e+87,
  1.268869321858842e+89,  8.247650592082472e+90,  5.443449390774431e+92,  3.647111091818868e+94,
  2.480035542436831e+96,  1.711224524281413e+98,  1.197857166996989e+100, 8.504785885678624e+101,
  6.123445837688610e+103, 4.470115461512685e+105, 3.307885441519387e+107, 2.480914081139540e+109,
  1.885494701666050e+111, 1.451830920282859e+113, 1.132428117820630e+115, 8.946182130782977e+116,
  7.156945704626381e+118, 5.797126020747368e+120, 4.753643337012842e+122, 3.945523969720659e+124,
  3.314240134565353e+126, 2.817104114380551e+128, 2.422709538367274e+130, 2.107757298379528e+132,
  1.854826422573984e+134, 1.650795516090846e+136, 1.485715964481762e+138, 1.352001527678403e+140,
  1.243841405464131e+142, 1.156772507081642e+144, 1.087366156656743e+146, 1.032997848823906e+148,
  9.916779348709498e+149, 9.619275968248212e+151, 9.426890448883246e+153, 9.332621544394417e+155,
  9.332621544394415e+157, 9.425947759838358e+159, 9.614466715035127e+161, 9.902900716486180e+163,
  1.029901674514563e+166, 1.081396758240291e+168, 1.146280563734709e+170, 1.226520203196138e+172,
  1.324641819451829e+174, 1.443859583202494e+176, 1.588245541522743e+178, 1.762952551090245e+180,
  1.974506857221074e+182, 2.231192748659814e+184, 2.543559733472188e+186, 2.925093693493017e+188,
  3.393108684451899e+190, 3.969937160808721e+192, 4.684525849754290e+194, 5.574585761207606e+196,
  6.689502913449128e+198, 8.094298525273445e+200, 9.875044200833603e+202, 1.214630436702533e+205,
  1.506141741511141e+207, 1.882677176888926e+209, 2.372173242880047e+211, 3.012660018457660e+213,
  3.856204823625805e+215, 4.974504222477287e+217, 6.466855489220474e+219, 8.471580690878822e+221,
  1.118248651196005e+224, 1.487270706090686e+226, 1.992942746161519e+228, 2.690472707318050e+230,
  3.659042881952549e+232, 5.012888748274991e+234, 6.917786472619488e+236, 9.615723196941092e+238,
  1.346201247571753e+241, 1.898143759076171e+243, 2.695364137888163e+245, 3.854370717180073e+247,
  5.550293832739306e+249, 8.047926057471993e+251, 1.174997204390911e+254, 1.727245890454639e+256,
  2.556323917872867e+258, 3.808922637630570e+260, 5.713383956445855e+262, 8.627209774233242e+264,
  1.311335885683453e+267, 2.006343905095682e+269, 3.089769613847352e+271, 4.789142901463394e+273,
  7.471062926282895e+275, 1.172956879426415e+278, 1.853271869493735e+280, 2.946702272495038e+282,
  4.714723635992062e+284, 7.590705053947219e+286, 1.229694218739449e+289, 2.004401576545302e+291,
  3.287218585534297e+293, 5.423910666131590e+295, 9.003691705778440e+297, 1.503616514864999e+300,
  2.526075744973198e+302, 4.269068009004706e+304, 7.257415615308000e+306};

//////////////////////////////// used in W3JS finish



// I/O routines

void save_vec(Vec &m, const char * filename, const char * name)
{
  std::cout << "Saving vector " << name << " in " << filename << "..." << std::endl;

  PetscErrorCode ierr;
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer); CHKERRV(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject)m,name); CHKERRV(ierr);
  ierr = VecView(m,viewer); CHKERRV(ierr);
  PetscViewerDestroy(&viewer);
}

void save_mat(Mat &m, const char * filename, const char * name)
{
  std::cout << "Saving matrix " << name << " in " << filename << "..." << std::endl;

  PetscErrorCode ierr;
  PetscViewer viewer;
  ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,filename,&viewer); CHKERRV(ierr);
  ierr = PetscViewerPushFormat(viewer,PETSC_VIEWER_ASCII_MATLAB); CHKERRV(ierr);
  ierr = PetscObjectSetName((PetscObject)m,name); CHKERRV(ierr);
  ierr = MatView(m,viewer); CHKERRV(ierr);
  PetscViewerDestroy(&viewer);
}

void read_vec(std::string filename, std::vector<double> &vec)
{
	int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	if (mpi_rank == 0) std::cout << "Reading vector data from " << filename << std::endl;

	std::ifstream myFile(filename);
	std::string line;	

  if (not myFile.good()) std::cerr << "\nERROR: File " << filename << " does not exist!\n" << std::endl;

	while(getline(myFile, line))
	{
    std::istringstream lineStream(line);
    double first;
    lineStream >> first;
    vec.push_back(first);		
	} 
}

void print_vec(const std::vector<double> &vec)
{
  for (auto x : vec) std::cout << x << std::endl;
}

void print_local_sizes(const Mat &M){
    
  PetscErrorCode ierr;
  PetscInt m, n;
  ierr = MatGetLocalSize(M, &m, &n);CHKERRV(ierr);
  std::cout << "m_local = " << m << std::endl;
  std::cout << "n_local = " << n << std::endl;
}

void print_global_sizes(const Mat &M){
    
  PetscErrorCode ierr;
  PetscInt m, n;
  ierr = MatGetSize(M, &m, &n);CHKERRV(ierr);
  std::cout << "m_global = " << m << std::endl;
  std::cout << "n_global = " << n << std::endl;
}

void create_identity_matrix(int size, Mat &Id)
{
  PetscErrorCode ierr;

  ierr = MatCreate(PETSC_COMM_WORLD, &Id);CHKERRV(ierr);  
  ierr = MatSetType(Id, MATAIJ);
  ierr = MatSetSizes(Id,PETSC_DECIDE,PETSC_DECIDE,size,size);CHKERRV(ierr);
  ierr = MatSeqAIJSetPreallocation(Id,1,NULL);CHKERRV(ierr);
  ierr = MatMPIAIJSetPreallocation(Id,1,NULL,0,NULL);CHKERRV(ierr);
  ierr = MatAssemblyBegin(Id,MAT_FINAL_ASSEMBLY);CHKERRV(ierr); 
  ierr = MatAssemblyEnd(Id,MAT_FINAL_ASSEMBLY);CHKERRV(ierr); 
  ierr = MatShift(Id, 1.0);CHKERRV(ierr);
}

std::vector<double> assemble_propagation_matrix(const std::vector<double> &etas_and_rhos)
{
  // check input sizes
  if (etas_and_rhos.size() != 7) std::cout << "\nERROR in assemble_propagation_matrix()\n" << std::endl;

  // etas_and_rhos = [eta0 eta1 eta2 eta3 rho1 rho2 rho3]

  std::vector<double> K(16);

  // etas
  K[0]  = etas_and_rhos[0];
  K[5]  = etas_and_rhos[0];
  K[10] = etas_and_rhos[0];
  K[15] = etas_and_rhos[0];
  K[1]  = etas_and_rhos[1];
  K[4]  = etas_and_rhos[1];
  K[2]  = etas_and_rhos[2];
  K[8]  = etas_and_rhos[2];
  K[3]  = etas_and_rhos[3];
  K[12] = etas_and_rhos[3];

  // rhos    
  K[6]  =  etas_and_rhos[6];
  K[9]  = -etas_and_rhos[6];
  K[7]  = -etas_and_rhos[5];
  K[13] =  etas_and_rhos[5];
  K[11] =  etas_and_rhos[4];
  K[14] = -etas_and_rhos[4];
      
  return K;
}

// entries divided by eta_I
std::vector<double> assemble_propagation_matrix_scaled(const std::vector<double> &etas_and_rhos)
{
  // check input sizes
  if (etas_and_rhos.size() != 7) std::cout << "\nERROR in assemble_propagation_matrix()\n" << std::endl;

  // etas_and_rhos = [eta0 eta1 eta2 eta3 rho1 rho2 rho3]

  std::vector<double> K(16);

  if ( etas_and_rhos[0] == 0) std::cout << "\n WARNING: eta_I = 0!\n" << std::endl;
  
  const double scaling_factor = 1.0 / etas_and_rhos[0];

  // etas
  K[0]  = 1.0;
  K[5]  = 1.0;
  K[10] = 1.0;
  K[15] = 1.0;
  K[1]  = scaling_factor * etas_and_rhos[1];
  K[4]  = scaling_factor * etas_and_rhos[1];
  K[2]  = scaling_factor * etas_and_rhos[2];
  K[8]  = scaling_factor * etas_and_rhos[2];
  K[3]  = scaling_factor * etas_and_rhos[3];
  K[12] = scaling_factor * etas_and_rhos[3];

  // rhos    
  K[6]  =  scaling_factor * etas_and_rhos[6];
  K[9]  = -scaling_factor * etas_and_rhos[6];
  K[7]  = -scaling_factor * etas_and_rhos[5];
  K[13] =  scaling_factor * etas_and_rhos[5];
  K[11] =  scaling_factor * etas_and_rhos[4];
  K[14] = -scaling_factor * etas_and_rhos[4];
          
  return K;
}


std::vector<double> assemble_propagation_matrix(const std::vector<double> &etas, const std::vector<double> &rhos){

  // check input sizes
  if (etas.size() != 4 or rhos.size() != 4) std::cout << "\nERROR in assemble_propagation_matrix()\n" << std::endl;
  
  if ( etas[0] == 0) std::cout << "\n WARNING: eta_I = 0!\n" << std::endl;

  std::vector<double> K(16);

  // etas
  K[0]  = etas[0];
  K[5]  = etas[0];
  K[10] = etas[0];
  K[15] = etas[0];
  K[1]  = etas[1];
  K[4]  = etas[1];
  K[2]  = etas[2];
  K[8]  = etas[2];
  K[3]  = etas[3];
  K[12] = etas[3];

  // rhos    
  K[6]  =  rhos[3];
  K[9]  = -rhos[3];
  K[7]  = -rhos[2];
  K[13] =  rhos[2];
  K[11] =  rhos[1];
  K[14] = -rhos[1];
      
  return K;

}

std::vector<double> assemble_propagation_matrix_scaled(const std::vector<double> &etas, const std::vector<double> &rhos)
{
  // check input sizes
  if (etas.size() != 4 or rhos.size() != 4) std::cout << "\nERROR in assemble_propagation_matrix()\n" << std::endl;

  if ( etas[0] == 0) std::cout << "\n WARNING: eta_I = 0!\n" << std::endl;

  std::vector<double> K(16);
  
  const double scaling_factor = 1.0 / etas[0];
  
  // etas
  K[0]  = 1.0;
  K[5]  = 1.0;
  K[10] = 1.0;
  K[15] = 1.0;
  K[1]  = scaling_factor * etas[1];
  K[4]  = scaling_factor * etas[1];
  K[2]  = scaling_factor * etas[2];
  K[8]  = scaling_factor * etas[2];
  K[3]  = scaling_factor * etas[3];
  K[12] = scaling_factor * etas[3];

  // rhos    
  K[6]  =  scaling_factor * rhos[3];
  K[9]  = -scaling_factor * rhos[3];
  K[7]  = -scaling_factor * rhos[2];
  K[13] =  scaling_factor * rhos[2];
  K[11] =  scaling_factor * rhos[1];
  K[14] = -scaling_factor * rhos[1];
  
  return K;
}

void print_propagation_matrix(const std::vector<double> &K)
{

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
  
  if (mpi_rank == 0)
  {      
    if (K.size() != 16) std::cout << "\nERROR in print_propagation_matrix()\n" << std::endl;

    for (int i = 0; i < 16; ++i)
    {
      if (i % 4 == 0) std::cout << std::endl;
      
      std::cout << K[i] << " ";
    }

    std::cout << std::endl;
  }
}

void print_Stokes(const std::vector<double> &I){

  int mpi_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

  if (mpi_rank == 0)
  {
    for (int i = 0; i < 4; ++i)
    {
      std::cout << I[i] << " ";
    }

    std::cout << std::endl;
  }
}

// W3JS from Simone 
double W3JS(int J1, int J2, int J3, int M1, int M2, int M3) {

  int IA, IB, IC, ID, IE, IF, IG, IH;
  int JSUM, Z, ZMIN, ZMAX;
  double CC, CC1, CC2, DENOM;

  // // correction from Simone
  // J1 *= 2;
  // J2 *= 2;
  // J3 *= 2;
  // M1 *= 2;
  // M2 *= 2;
  // M3 *= 2;
  // /////////////////////////
  

  if (M1 + M2 + M3 != 0)
    return 0.0;
  IA = J1 + J2;
  if (J3 > IA)
    return 0.0;  // disuguagliaza triangilare prima riga.
  IB = J1 - J2;
  if (J3 < abs(IB))
    return 0.0;  // minore della differenza ritorna 0
  JSUM = J3 + IA;
  IC   = J1 - M1;
  ID   = J2 - M2;
  if (abs(M1) > J1 || abs(M2) > J2 || abs(M3) > J3)
    return 0.0;

  IE   = J3 - J2 + M1;
  IF   = J3 - J1 - M2;
  ZMIN = MAX3(0, -IE, -IF);
  IG   = IA - J3;
  IH   = J2 + M2;
  ZMAX = MIN3(IG, IH, IC);
  CC   = 0.0;

  for (Z = ZMIN; Z <= ZMAX; Z += 2) {
    DENOM = g_fact[Z / 2] * g_fact[(IG - Z) / 2] * g_fact[(IC - Z) / 2] * g_fact[(IH - Z) / 2] * g_fact[(IE + Z) / 2] *
            g_fact[(IF + Z) / 2];
    if (Z % 4)
      DENOM = -DENOM;
    CC = CC + 1.0 / DENOM;
  }

  CC1 = g_fact[IG / 2] * g_fact[(J3 + IB) / 2] * g_fact[(J3 - IB) / 2] / g_fact[(JSUM + 2) / 2];
  CC2 = g_fact[(J1 + M1) / 2] * g_fact[IC / 2] * g_fact[IH / 2] * g_fact[ID / 2] * g_fact[(J3 - M3) / 2] *
        g_fact[(J3 + M3) / 2];
  CC *= sqrt(CC1 * CC2);
  if (Mod(IB - M3, 4))
    CC = -CC;

  if (fabs(CC) < 1e-8)
    return 0;
  else
    return CC;
}

 
std::vector<double> solve_4_by_4_system(const std::vector<double>  &K, const std::vector<double> &rhs)
{
  // check input sizes
  if (K.size() != 16 || rhs.size() != 4 ) std::cout << "\nERROR in solve_4_by_4_system()\n" << std::endl;
              
  std::vector<double> sol(4);
  std::vector<double> inv(16);

  inv[0] = K[5]  * K[10] * K[15] - 
           K[5]  * K[11] * K[14] - 
           K[9]  * K[6]  * K[15] + 
           K[9]  * K[7]  * K[14] +
           K[13] * K[6]  * K[11] - 
           K[13] * K[7]  * K[10];

  inv[4] = -K[4]  * K[10] * K[15] + 
            K[4]  * K[11] * K[14] + 
            K[8]  * K[6]  * K[15] - 
            K[8]  * K[7]  * K[14] - 
            K[12] * K[6]  * K[11] + 
            K[12] * K[7]  * K[10];

  inv[8] = K[4]  * K[9] * K[15] - 
           K[4]  * K[11] * K[13] - 
           K[8]  * K[5] * K[15] + 
           K[8]  * K[7] * K[13] + 
           K[12] * K[5] * K[11] - 
           K[12] * K[7] * K[9];

  inv[12] = -K[4]  * K[9] * K[14] + 
             K[4]  * K[10] * K[13] +
             K[8]  * K[5] * K[14] - 
             K[8]  * K[6] * K[13] - 
             K[12] * K[5] * K[10] + 
             K[12] * K[6] * K[9];

  inv[1] = -K[1]  * K[10] * K[15] + 
            K[1]  * K[11] * K[14] + 
            K[9]  * K[2] * K[15] - 
            K[9]  * K[3] * K[14] - 
            K[13] * K[2] * K[11] + 
            K[13] * K[3] * K[10];

  inv[5] = K[0]  * K[10] * K[15] - 
           K[0]  * K[11] * K[14] - 
           K[8]  * K[2] * K[15] + 
           K[8]  * K[3] * K[14] + 
           K[12] * K[2] * K[11] - 
           K[12] * K[3] * K[10];

  inv[9] = -K[0]  * K[9] * K[15] + 
            K[0]  * K[11] * K[13] + 
            K[8]  * K[1] * K[15] - 
            K[8]  * K[3] * K[13] - 
            K[12] * K[1] * K[11] + 
            K[12] * K[3] * K[9];

  inv[13] = K[0]  * K[9] * K[14] - 
            K[0]  * K[10] * K[13] - 
            K[8]  * K[1] * K[14] + 
            K[8]  * K[2] * K[13] + 
            K[12] * K[1] * K[10] - 
            K[12] * K[2] * K[9];

  inv[2] = K[1]  * K[6] * K[15] - 
           K[1]  * K[7] * K[14] - 
           K[5]  * K[2] * K[15] + 
           K[5]  * K[3] * K[14] + 
           K[13] * K[2] * K[7] - 
           K[13] * K[3] * K[6];

  inv[6] = -K[0]  * K[6] * K[15] + 
            K[0]  * K[7] * K[14] + 
            K[4]  * K[2] * K[15] - 
            K[4]  * K[3] * K[14] - 
            K[12] * K[2] * K[7] + 
            K[12] * K[3] * K[6];

  inv[10] = K[0]  * K[5] * K[15] - 
            K[0]  * K[7] * K[13] - 
            K[4]  * K[1] * K[15] + 
            K[4]  * K[3] * K[13] + 
            K[12] * K[1] * K[7] - 
            K[12] * K[3] * K[5];

  inv[14] = -K[0]  * K[5] * K[14] + 
             K[0]  * K[6] * K[13] + 
             K[4]  * K[1] * K[14] - 
             K[4]  * K[2] * K[13] - 
             K[12] * K[1] * K[6] + 
             K[12] * K[2] * K[5];

  inv[3] = -K[1] * K[6] * K[11] + 
            K[1] * K[7] * K[10] + 
            K[5] * K[2] * K[11] - 
            K[5] * K[3] * K[10] - 
            K[9] * K[2] * K[7] + 
            K[9] * K[3] * K[6];

  inv[7] = K[0] * K[6] * K[11] - 
           K[0] * K[7] * K[10] - 
           K[4] * K[2] * K[11] + 
           K[4] * K[3] * K[10] + 
           K[8] * K[2] * K[7] - 
           K[8] * K[3] * K[6];

  inv[11] = -K[0] * K[5] * K[11] + 
             K[0] * K[7] * K[9] + 
             K[4] * K[1] * K[11] - 
             K[4] * K[3] * K[9] - 
             K[8] * K[1] * K[7] + 
             K[8] * K[3] * K[5];

  inv[15] = K[0] * K[5] * K[10] - 
            K[0] * K[6] * K[9] - 
            K[4] * K[1] * K[10] + 
            K[4] * K[2] * K[9] + 
            K[8] * K[1] * K[6] - 
            K[8] * K[2] * K[5];

  double det = K[0] * inv[0] + K[1] * inv[4] + K[2] * inv[8] + K[3] * inv[12];

  if (det == 0) std::cout << "\nERROR: K matrix not invertible!\n" << std::endl;    

  // if (std::abs(det) < 1e-5) std::cout << "\nWARNING: K matrix is ill conditioned! det(K) = " << std::abs(det)  << "\n" << std::endl;    

  det = 1.0 / det;

  double sum;

  // compute inv(A) * rhs
  for (int i = 0; i < 4; ++i)
  {
      sum = 0;

      for (int j = 0; j < 4; ++j)
      {
          sum += inv[4 * i + j] * rhs[j];
      }

      sol[i] = det * sum;
  }
  
  return sol;
}


// for matrix K having identity diagonal
std::vector<double> solve_4_by_4_system_optimized(const std::vector<double>  &K, const std::vector<double> &rhs)
{
  // check input sizes
  if (K.size() != 16 || rhs.size() != 4 ) std::cout << "\nERROR in solve_4_by_4_system()\n" << std::endl;
              
  std::vector<double> sol(4);
  std::vector<double> inv(16);

  if (K[0] != 1.0 or K[5] != 1.0 or K[10] != 1.0 or K[15] != 1.0)  std::cout << "\nERROR can not use solve_4_by_4_system_optimized()!\n" << std::endl;

  inv[0] = 1.0 - 
           K[11] * K[14] - 
           K[9]  * K[6]  + 
           K[9]  * K[7]  * K[14] +
           K[13] * K[6]  * K[11] - 
           K[13] * K[7] ;

  inv[4] = -K[4] + 
            K[4]  * K[11] * K[14] + 
            K[8]  * K[6]  - 
            K[8]  * K[7]  * K[14] - 
            K[12] * K[6]  * K[11] + 
            K[12] * K[7];

  inv[8] = K[4]  * K[9] - 
           K[4]  * K[11] * K[13] - 
           K[8]  + 
           K[8]  * K[7] * K[13] + 
           K[12] * K[11] - 
           K[12] * K[7] * K[9];

  inv[12] = -K[4]  * K[9] * K[14] + 
             K[4]  * K[13] +
             K[8]  * K[14] - 
             K[8]  * K[6] * K[13] - 
             K[12] + 
             K[12] * K[6] * K[9];

  inv[1] = -K[1] + 
            K[1]  * K[11] * K[14] + 
            K[9]  * K[2] - 
            K[9]  * K[3] * K[14] - 
            K[13] * K[2] * K[11] + 
            K[13] * K[3];

  inv[5] = 1.0 - 
           K[11] * K[14] - 
           K[8]  * K[2] + 
           K[8]  * K[3] * K[14] + 
           K[12] * K[2] * K[11] - 
           K[12] * K[3];

  inv[9] = - K[9] + 
            K[11] * K[13] + 
            K[8]  * K[1] - 
            K[8]  * K[3] * K[13] - 
            K[12] * K[1] * K[11] + 
            K[12] * K[3] * K[9];

  inv[13] = K[9] * K[14] - 
            K[13] - 
            K[8]  * K[1] * K[14] + 
            K[8]  * K[2] * K[13] + 
            K[12] * K[1] - 
            K[12] * K[2] * K[9];

  inv[2] = K[1] * K[6] - 
           K[1] * K[7] * K[14] - 
           K[2] + 
           K[3] * K[14] + 
           K[13] * K[2] * K[7] - 
           K[13] * K[3] * K[6];

  inv[6] = -K[6] + 
            K[7] * K[14] + 
            K[4]  * K[2] - 
            K[4]  * K[3] * K[14] - 
            K[12] * K[2] * K[7] + 
            K[12] * K[3] * K[6];

  inv[10] = 1.0 - 
            K[7] * K[13] - 
            K[4]  * K[1] + 
            K[4]  * K[3] * K[13] + 
            K[12] * K[1] * K[7] - 
            K[12] * K[3];

  inv[14] = - K[14] + 
             K[6] * K[13] + 
             K[4]  * K[1] * K[14] - 
             K[4]  * K[2] * K[13] - 
             K[12] * K[1] * K[6] + 
             K[12] * K[2];

  inv[3] = -K[1] * K[6] * K[11] + 
            K[1] * K[7] + 
            K[2] * K[11] - 
            K[3] - 
            K[9] * K[2] * K[7] + 
            K[9] * K[3] * K[6];

  inv[7] = K[6] * K[11] - 
           K[7]- 
           K[4] * K[2] * K[11] + 
           K[4] * K[3] + 
           K[8] * K[2] * K[7] - 
           K[8] * K[3] * K[6];

  inv[11] = -K[11] + 
             K[7] * K[9] + 
             K[4] * K[1] * K[11] - 
             K[4] * K[3] * K[9] - 
             K[8] * K[1] * K[7] + 
             K[8] * K[3];

  inv[15] = 1.0 - 
            K[6] * K[9] - 
            K[4] * K[1] + 
            K[4] * K[2] * K[9] + 
            K[8] * K[1] * K[6] - 
            K[8] * K[2] ;

  double det = inv[0] + K[1] * inv[4] + K[2] * inv[8] + K[3] * inv[12];

  if (det == 0) std::cout << "\nERROR: K matrix not invertible!\n" << std::endl;    

  // if (std::abs(det) < 1e-5) std::cout << "\nWARNING: K matrix is ill conditioned! det(K) = " << std::abs(det)  << "\n" << std::endl;    

  det = 1.0 / det;

  double sum;

  // compute inv(A) * rhs
  for (int i = 0; i < 4; ++i)
  {
      sum = 0;

      for (int j = 0; j < 4; ++j)
      {
          sum += inv[4 * i + j] * rhs[j];
      }

      sol[i] = det * sum;
  }
  
  return sol;
}


std::vector<double> refine_vector(const std::vector<double> &v)
{
  const size_t n = v.size();
  const size_t n_fn = 2 * n - 1;

  std::vector<double> v_fine;

  v_fine.reserve(n_fn);
  v_fine.resize(n_fn);  

  // copy
  for (size_t i = 0; i < n_fn; i = i + 2)
  {
    v_fine[i] = v[i/2];
  }

  // fill new values
  for (size_t i = 1; i < n_fn - 1; i = i + 2)
  {
    v_fine[i] = 0.5 * (v[i/2] + v[i/2 + 1]);
  }

  return v_fine;
}


std::vector<double> refine_vector_blocked(const std::vector<double> &v, const size_t block_size)
{

  const size_t n = v.size();
  
  const size_t number_of_blocks = n / block_size;

  if (n % block_size > 0) std::cout << "WARNING: error in block_size!" << std::endl;
  
  const size_t n_fn = (2 * number_of_blocks - 1 ) * block_size;

  std::vector<double> v_fine;

  v_fine.reserve(n_fn);
  v_fine.resize(n_fn);  

  // copy
  for (size_t i = 0; i < n_fn; i = i + 2 * block_size)
  {
    for (size_t j = 0; j < block_size; ++j)
    {
      v_fine[i + j] = v[i/2 + j]; 
    }
  }

  // fill
  size_t cr_index;
  size_t i_over_bs;

  for (size_t i = block_size; i < n_fn - block_size; i = i + 2 * block_size)
  {
    i_over_bs = i/block_size;

    for (size_t j = 0; j < block_size; ++j)
    {
      cr_index = block_size * (i_over_bs/2) + j;

      v_fine[i + j] = 0.5 * (v[cr_index] + v[cr_index + block_size]);
    }    
  }

  return v_fine;
}

// refine lines 
std::vector<double> refine_vector_blocked2(const std::vector<double> &v, const size_t block_size)
{
  const size_t n = v.size();
  const size_t number_of_blocks = n / block_size;
  const size_t block_size_fn = 2 * block_size - 1;
  const size_t n_fn = number_of_blocks * block_size_fn;

  
  std::vector<double> v_fine;

  v_fine.reserve(n_fn);
  v_fine.resize(n_fn);  

  size_t index;

  for (int i = 0; i < (int)number_of_blocks; ++i)
  {

    // copy
    for (size_t j = 0; j < block_size_fn; j = j + 2)
    {

      index = i * block_size_fn + j;

      v_fine[index] = v[index/2];
    }

    // fill new values
    for (size_t j = 1; j < block_size_fn - 1; j = j + 2)
    {

      index = i * block_size_fn + j;

      v_fine[index] = 0.5 * (v[index/2] + v[index/2 + 1]);
    }
    
  }

  return v_fine;
}


double pow_gen(const double x, const double exponent){

  if (x >= 0)
  {
    return std::pow(x,exponent);
  }
  else
  {
    return -std::pow(-x,exponent);
  }
}








#include "bp.hpp"

#include <thread>
#include "linalg/simd.hpp"
using Real = AutoSIMDTraits<double,double>::vreal_t;
#define SIMD_SIZE (MFEM_SIMD_BYTES/sizeof(double))

template<int DIM, int DX0, int DX1> inline static
void KSetup1(const int ndofs,
		const int vdim, const int NE,
		const double * __restrict__ J0,
		const double * __restrict__ w,
		double * __restrict__ dx) {
	assert(vdim == 1);
	static constexpr int Q1D = 5;
	dbg("SIMD_SIZE:%d", SIMD_SIZE);

	// kernel operations: u,G,*,D,*,v,G,T,*,Ye

	const auto J = Reshape(J0, DIM,DIM, Q1D,Q1D,Q1D, NE);
	const auto W = Reshape(w, Q1D,Q1D,Q1D);
	auto DX = Reshape(dx, DX0,DX1, Q1D,Q1D,Q1D, NE);

	MFEM_FORALL_3D(e, NE, Q1D, Q1D, Q1D,{
	MFEM_FOREACH_THREAD(qz,z,Q1D){
		MFEM_FOREACH_THREAD(qy,y,Q1D){
			MFEM_FOREACH_THREAD(qx,x,Q1D){
				const double irw = W(qx,qy,qz);
				const double *Jtr = &J(0,0,qx,qy,qz,e);
				const double detJ = kernels::Det<DIM>(Jtr);
				const double wd = irw * detJ;
				double Jrt[DIM*DIM];
				kernels::CalcInverse<DIM>(Jtr, Jrt);
				double A[DX0*DX1];
				double D[DX0*DX1] = {wd,0,0,0,wd,0,0,0,wd};
				kernels::MultABt(DIM,DIM,DIM,D,Jrt,A);
				kernels::Mult(DIM,DIM,DIM,A,Jrt,&DX(0,0,qx,qy,qz,e));
			}
		}
	}
		MFEM_SYNC_THREAD;
	});
}

template<int DIM, int DX0, int DX1> inline static
void KMult1(const int ndofs, const int vdim, const int NE,
		const double * __restrict__ B,
		const double * __restrict__ G,
		const int * __restrict__ map,
		const double * __restrict__ dx,
		const double * __restrict__ xd,
		double * __restrict__ yd) {

	static constexpr int D1D = 4;
	static constexpr int MD1 = 4;
	static constexpr int Q1D = 5;
	static constexpr int MQ1 = 5;
	static constexpr int SMS = SIMD_SIZE;

	// kernel operations: u,G,*,D,*,v,G,T,*,Ye

	assert(vdim == 1);

	const auto b = Reshape(B, Q1D, D1D);
	const auto g = Reshape(G, Q1D, D1D);
	const auto DX = Reshape(dx, DX0,DX1, Q1D,Q1D,Q1D, NE);
	const auto MAP = Reshape(map, D1D,D1D,D1D, NE);
	const auto XD = Reshape(xd, ndofs/*, vdim*/);
	auto YD = Reshape(yd, ndofs/*, vdim*/);

	double BG[2][MQ1*MD1];
	kernels::LoadBG<MD1,MQ1>(D1D,Q1D,b,g,BG);

	double BGt[2][MQ1*MD1];
	kernels::LoadBGt<MD1,MQ1>(D1D,Q1D,b,g,BGt);

#ifndef MFEM_USE_THREADS
	//#pragma omp parallel for
	MFEM_VERIFY((NE % SIMD_SIZE) == 0, "NE vs SIMD_SIZE error!")
#define BATCH_SIZE 32
	MFEM_VERIFY((NE % BATCH_SIZE) == 0, "NE vs BATCH_SIZE error!")
	MFEM_VERIFY(((NE/BATCH_SIZE) % SIMD_SIZE) == 0, "NE/BATCH_SIZE vs SIMD_SIZE error!")
	for(size_t eb = 0; eb < (NE/(BATCH_SIZE*SIMD_SIZE)); eb+=1) {
	for(size_t e = eb*BATCH_SIZE*SIMD_SIZE; e < (eb+1)*BATCH_SIZE*SIMD_SIZE; e+=SIMD_SIZE) {
#else
	std::vector<std::thread> threads;
	static const unsigned int num_threads = std::thread::hardware_concurrency();
	dbg("NE:%d, num_threads:%d",NE, num_threads);
	MFEM_VERIFY((NE % num_threads) == 0, "NE vs #Threads error")
	int e0 = 0;
	const int NEB = NE / num_threads;
	int NE0 = NEB;
	for(unsigned int tid = 0; tid < num_threads; ++tid) {
		threads.push_back(std::thread(
	[&](const int tid, const int e0, const int NE0){
	//printf("[#%d] e0:%d, NE0:%d", tid, e0, NE0);
	for(size_t e = e0; e < NE0; e+=SIMD_SIZE) {
#endif
		Real DDD[MD1*MD1*MD1];
		kernels::LoadXDGather<MD1,Real,SMS>(e,D1D,MAP,XD,DDD);
		// kernel operations: u,G,*,D,*,v,G,T,*,Ye
		// [push] trial u:2
		Real sm0[3][MQ1*MQ1*MQ1];
		Real sm1[3][MQ1*MQ1*MQ1];
		Real (*DDQ)[MD1*MD1*MQ1] = (Real (*)[MD1*MD1*MQ1]) (sm0);
		Real (*DQQ)[MD1*MQ1*MQ1] = (Real (*)[MD1*MQ1*MQ1]) (sm1);
		Real (*QQQ)[MQ1*MQ1*MQ1] = (Real (*)[MQ1*MQ1*MQ1]) (sm0);
		// Grad(u)
		kernels::Grad1X<MD1,MQ1>(D1D,Q1D,BG,DDD,DDQ);
		kernels::Grad1Y<MD1,MQ1>(D1D,Q1D,BG,DDQ,DQQ);
		kernels::Grad1Z<MD1,MQ1>(D1D,Q1D,BG,DQQ,QQQ);
		// [ pop] u
	for(int qz = 0; qz < Q1D; qz++){
		for(int qy = 0; qy < Q1D; qy++){
			for(int qx = 0; qx < Q1D; qx++){
				Real u[DX0], v[DX0];
				kernels::PullGrad1<MQ1>(qx,qy,qz,QQQ,u);
				kernels::Mult(DX0,DX1,&DX(0,0,qx,qy,qz,e),u,v);
				kernels::PushGrad1<MQ1>(qx,qy,qz,v,QQQ);
		// [push] test v:2
			}
		}
	}
		// Grad(v)
		kernels::Grad1Zt<MD1,MQ1>(D1D,Q1D,BGt,QQQ,DQQ);
		kernels::Grad1Yt<MD1,MQ1>(D1D,Q1D,BGt,DQQ,DDQ);
		kernels::Grad1XtDScatter<MD1,MQ1,Real,SMS>(D1D,Q1D,BGt,DDQ,MAP,YD,e);
		// [ pop] v
	}} // Element for loop
#ifdef MFEM_USE_THREADS
	}, tid, e0, NE0)); // lambda & thread vector push_back 
	e0 += NEB; NE0 += NEB;
	} // Thread for loop
	for(auto &thr : threads) { thr.join(); }
#endif
} // KMult1

#ifndef GEOM
#define GEOM Geometry::CUBE
#endif

#ifndef MESH_P
#define MESH_P 1
#endif

#ifndef SOL_P
#define SOL_P 3
#endif

#ifndef IR_ORDER
#define IR_ORDER 0
#endif

#ifndef IR_TYPE
#define IR_TYPE 0
#endif

#ifndef PROBLEM
#define PROBLEM 0
#endif

#ifndef VDIM
#define VDIM 1
#endif

int main(int argc, char* argv[]){
	int status = 0;
	int num_procs, myid;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	assert(VDIM==1); assert(MESH_P==1); assert(IR_TYPE==0); assert(IR_ORDER==0); assert(PROBLEM==0); assert(GEOM==Geometry::CUBE); const char *mesh_file = "../../data/hex-01x01x01.mesh"; int ser_ref_levels = 3; int par_ref_levels = 1; Array<int> nxyz; int order = SOL_P; const char *basis_type = "G"; bool static_cond = false; const char *pc = "none"; bool perf = true; bool matrix_free = true; int max_iter = 50; bool visualization = 0; OptionsParser args(argc, argv); args.AddOption(&mesh_file, "-m", "--mesh", "Mesh file to use."); args.AddOption(&ser_ref_levels, "-rs", "--refine-serial", "Number of times to refine the mesh uniformly in serial."); args.AddOption(&par_ref_levels, "-rp", "--refine-parallel", "Number of times to refine the mesh uniformly in parallel."); args.AddOption(&nxyz, "-c", "--cartesian-partitioning", "Use Cartesian partitioning."); args.AddOption(&order, "-o", "--order", "Finite element order (polynomial degree) or -1 for" " isoparametric space."); args.AddOption(&basis_type, "-b", "--basis-type", "Basis: G - Gauss-Lobatto, P - Positive, U - Uniform"); args.AddOption(&perf, "-perf", "--hpc-version", "-std", "--standard-version", "Enable high-performance, tensor-based, assembly/evaluation."); args.AddOption(&matrix_free, "-mf", "--matrix-free", "-asm", "--assembly", "Use matrix-free evaluation or efficient matrix assembly in " "the high-performance version."); args.AddOption(&pc, "-pc", "--preconditioner", "Preconditioner: lor - low-order-refined (matrix-free) AMG, " "ho - high-order (assembled) AMG, none."); args.AddOption(&static_cond, "-sc", "--static-condensation", "-no-sc", "--no-static-condensation", "Enable static condensation."); args.AddOption(&max_iter, "-mi", "--max-iter", "Maximum number of iterations."); args.AddOption(&visualization, "-vis", "--visualization", "-no-vis", "--no-visualization", "Enable or disable GLVis visualization."); args.Parse(); if (!args.Good()) { if (myid == 0) { args.PrintUsage(std::cout); } return 1; } if (myid == 0) { args.PrintOptions(std::cout); } ParMesh *pmesh = nullptr; { Mesh *mesh = new Mesh(mesh_file, 1, 1); int dim = mesh->Dimension(); { int ref_levels = (int)floor(log(10000./mesh->GetNE())/log(2.)/dim); ref_levels = (ser_ref_levels != -1) ? ser_ref_levels : ref_levels; for (int l = 0; l < ref_levels; l++) { if (myid == 0) { std::cout << "Serial refinement: level " << l << " -> level " << l+1 << " ..." << std::flush; } mesh->UniformRefinement(); MPI_Barrier(MPI_COMM_WORLD); if (myid == 0) { std::cout << " done." << std::endl; } } } MFEM_VERIFY(nxyz.Size() == 0 || nxyz.Size() == mesh->SpaceDimension(), "Expected " << mesh->SpaceDimension() << " integers with the " "option --cartesian-partitioning."); int *partitioning = nxyz.Size() ? mesh->CartesianPartitioning(nxyz) : NULL; pmesh = new ParMesh(MPI_COMM_WORLD, *mesh, partitioning); delete [] partitioning; delete mesh; { for (int l = 0; l < par_ref_levels; l++) { if (myid == 0) { std::cout << "Parallel refinement: level " << l << " -> level " << l+1 << " ..." << std::flush; } pmesh->UniformRefinement(); MPI_Barrier(MPI_COMM_WORLD); if (myid == 0) { std::cout << " done." << std::endl; } } } pmesh->PrintInfo(std::cout); }
	const int p = 3;
	const int dim = 3;
	auto &mesh = xfl::Mesh(pmesh);
	const int el = (dim==2)?xfl::quadrilateral:xfl::hexahedron;
	FiniteElementCollection *fe = xfl::FiniteElement("Lagrange", el, p);
	mfem::ParFiniteElementSpace *fes = xfl::FunctionSpace(mesh, fe);
	const Array<int> bc = xfl::DirichletBC(fes);
	xfl::Function x = xfl::Function(fes);
	x = 0.0;
	xfl::TrialFunction u = xfl::TrialFunction(fes);
	xfl::TestFunction v = xfl::TestFunction(fes);
	auto b = [&]() {
		constexpr const char *qs0 = "v";
		// var:[v], ops:[Xe,v,]
		// Test FES: 'fes':v (Eval)
		ParFiniteElementSpace *fes0 = fes;
		mfem::Operator *QM0 = nullptr;
		xfl::QForm QForm0(fes0, qs0, QM0);
		return QForm0;
	}();
	auto a = [&]() {
		constexpr const char *qs1 = "dot(grad(u), grad(v))";
		// var:[u,v], ops:[Xe,u,G,*,D,*,v,G,T,*,Ye]
		// Trial FES: 'fes':u (Grad)
		// Test FES: 'fes':v (Grad)
		ParFiniteElementSpace *fes1 = fes;
		struct QMult1: public xfl::Operator<3>{
			QMult1(const ParFiniteElementSpace *fes): xfl::Operator<3>(fes) { dbg(); Setup(); }
			~QMult1() { dbg(); }
			void Setup() {
				//dbg();
				dx.SetSize(NQ*NE*3*3, Device::GetDeviceMemoryType()); // DX shape: 3x3
				KSetup1<3,3,3>(NDOFS, VDIM, NE, J0.Read(), ir.GetWeights().Read(), dx.Write());
			}
			void Mult(const mfem::Vector &x, mfem::Vector &y) const {
				y = 0.0;
				KMult1<3,3,3>(NDOFS /*0*/,VDIM /*0*/, NE /*0*/,maps->B.Read(), maps->G.Read(), ER.GatherMap().Read(), dx.Read(), x.Read(), y.ReadWrite());
			}
		}; // QMult struct
		QMult1 *QM1 = new QMult1(fes);
		xfl::QForm QForm1(fes1, qs1, QM1);
		return QForm1;
	}();
	status |= xfl::benchmark(a==b, x, bc, 0, 50, -1);
	const bool glvis = true;
	if(glvis)xfl::plot(x);
	MPI_Finalize();
	return status;
}


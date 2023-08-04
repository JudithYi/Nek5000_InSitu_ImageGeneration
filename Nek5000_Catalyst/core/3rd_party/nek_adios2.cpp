#include "nek_adios2.h"

int lx1, ly1, lz1, lx2, ly2, lz2;
int lxyz1, lxyz2; 
int nelv, nelt;
int nelgv, nelgt;
int nelbv, nelbt;

std::size_t total1, start1, count1;
std::size_t total2, start2, count2;
std::size_t total3, start3, count3;
std::size_t total4, start4, count4;
std::size_t init_total, init_start, init_count;
adios2::ADIOS adios;
adios2::IO io;
std::vector<adios2::Engine> writer;
adios2::Variable<int> init_int_const;
adios2::Variable<double> init_double_const;
std::vector<int> vINT;
std::vector<double> vDOUBLE;
std::vector<int> vInSituCounter;
adios2::Variable<double> x;
adios2::Variable<double> y;
adios2::Variable<double> z;
adios2::Variable<double> p;
adios2::Variable<double> vx;
adios2::Variable<double> vy;
adios2::Variable<double> vz;
adios2::Variable<double> t;
adios2::Variable<double> bm1;
adios2::Variable<int> InsituCounter;
double dataTime = 0.0;
std::clock_t startT;
std::clock_t startTotal;
int rank, size, i;
bool if3d;
bool firstStep;
int writer_num = 0;
int insituCount = 1;
int world_rank, world_size, insitu_size;

void adios_writer_init(
    MPI_Comm & comm,
    MPI_Comm & worldComm, 
    std::string engineName
){  
    if(writer_num == 0){
	startTotal = std::clock();
        std::string configFile="config/config.xml";
        MPI_Comm_rank(worldComm, &world_rank);
        MPI_Comm_size(worldComm, &world_size);
        //MPI_Comm_split(worldComm, colour, world_rank, &comm);
        adios = adios2::ADIOS(configFile, comm);
	io = adios.DeclareIO("writer");
        int	worldComm_f = MPI_Comm_c2f(worldComm);
        int comm_f = MPI_Comm_c2f(comm);
        std::cout << "From writer[" << writer_num << "]: " << comm_f << " " << worldComm_f << " " << engineName << std::endl;
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
        insitu_size = world_size - size;
    }
    writer.resize(writer_num+1);
    writer[writer_num] = io.Open(engineName, adios2::Mode::Write,comm);
    //if(io.EngineType()=="InSituMPI")
    io.set_worldComm(&writer[writer_num],&worldComm);
}

extern "C" void adios2_setup_async_(
    const int *lx1_in,
    const int *ly1_in,
    const int *lz1_in,
    const int *lx2_in,
    const int *ly2_in,
    const int *lz2_in,
    const int *nelv_in,
    const int *nelt_in,
    const int *nelbv_in,
    const int *nelbt_in,
    const int *nelgv_in,
    const int *nelgt_in,
    const int *iostep_in,
    const int *if3d_in,
    const double *t_start_in,
    const double *dt_in,
    const double *xml_in,
    const double *yml_in,
    const double *zml_in,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *tem,
    const double *bm,
    const int *comm_int,
    const int *world_comm_int
){
    //writer0.BeginStep();
    if(writer_num == 0){
        lx1 = *lx1_in;
        ly1 = *ly1_in;
        lz1 = *lz1_in;
        lx2 = *lx2_in;
        ly2 = *ly2_in;
        lz2 = *lz2_in;
        nelgt = *nelgt_in;
        nelgv = *nelgv_in;
        if3d = true;
        if(*if3d_in==0) if3d = false;
        nelt = *nelt_in;
        nelv = *nelv_in;
        /* To make sure the in-situ core with rank i recieves the data from the simulation core with rank i*(size/insitu_size) to i*(size/insitu_size+1)-1 */ 
        //nelbv = *nelbv_in;
        //nelbt = *nelbt_in;
        MPI_Comm comm = MPI_Comm_f2c(*comm_int);
        
        std::vector<int>temp(size);
        std::vector<int>temp1(size);
        MPI_Allgather(&nelv, 1, MPI_INT, temp.data(), 1,MPI_INT, comm);
        MPI_Allgather(&nelt, 1, MPI_INT, temp1.data(), 1,MPI_INT, comm);
        nelbt = 0;
        nelbv = 0;
        for(i=0;i<rank;++i){
            nelbv += temp[i];
            nelbt += temp1[i];
        }
        
        init_count = 12;
        init_total = init_count*static_cast<std::size_t>(size);
        init_start = init_count*static_cast<std::size_t>(rank);
        vINT.resize(init_count);
        vINT[0] = lx1;
        vINT[1] = ly1;
        vINT[2] = lz1;
        vINT[3] = lx2;
        vINT[4] = ly2;
        vINT[5] = lz2;
        vINT[6] = nelv;
        vINT[7] = nelt;
        vINT[8] = nelgv;
        vINT[9] = nelgt;
        vINT[10] = *iostep_in;
        vINT[11] = *if3d_in;
        init_int_const = io.DefineVariable<int>("INT_CONST", {init_total}, {init_start}, {init_count});
        
        init_count = 2;
        init_total = init_count*static_cast<std::size_t>(size);
        init_start = init_count*static_cast<std::size_t>(rank);
        vDOUBLE.resize(init_count);
        vDOUBLE[0] = *t_start_in;
        vDOUBLE[1] = *dt_in;
        vInSituCounter.resize(1);

        init_double_const = io.DefineVariable<double>("DOUBLE_CONST", {init_total}, {init_start}, {init_count});
        lxyz1 = lx1 * ly1 * lz1;
        lxyz2 = lx2 * ly2 * lz2;
        total1 = static_cast<std::size_t>(lxyz1 * nelgv);
        start1 = static_cast<std::size_t>(lxyz1 * nelbv);
        count1 = static_cast<std::size_t>(lxyz1 * nelv);
        total2 = static_cast<std::size_t>(lxyz2 * nelgv);
        start2 = static_cast<std::size_t>(lxyz2 * nelbv);
        count2 = static_cast<std::size_t>(lxyz2 * nelv);
        total3 = static_cast<std::size_t>(lxyz1 * nelgt);
        start3 = static_cast<std::size_t>(lxyz1 * nelbt);
        count3 = static_cast<std::size_t>(lxyz1 * nelt);

        total4 = static_cast<std::size_t>(size);
        start4 = static_cast<std::size_t>(rank);
        count4 = 1;

        x = io.DefineVariable<double>("X", {total3}, {start3}, {count3});
        y = io.DefineVariable<double>("Y", {total3}, {start3}, {count3});
        z = io.DefineVariable<double>("Z", {total3}, {start3}, {count3});
        p = io.DefineVariable<double>("P", {total2}, {start2}, {count2});
        vx = io.DefineVariable<double>("VX", {total1}, {start1}, {count1});
        vy = io.DefineVariable<double>("VY", {total1}, {start1}, {count1});
        vz = io.DefineVariable<double>("VZ", {total1}, {start1}, {count1});
        t = io.DefineVariable<double>("T", {total3}, {start3}, {count3});
        bm1 = io.DefineVariable<double>("BM1", {total3}, {start3}, {count3});
        InsituCounter = io.DefineVariable<int>("InsituCount", {total4}, {start4}, {count4});
    }
    vInSituCounter[0] = insituCount;

    adios2::StepStatus status = writer[writer_num].BeginStep();
    if (status != adios2::StepStatus::OK){
	    std::cout << "Writer has not done begin step." << std::endl;
        writer[writer_num].EndStep();
        writer[writer_num].Close();
        //MPI_Finalize();
        return;
    }
    //writer = io.Open("globalArray", adios2::Mode::Write);
    //writer.BeginStep();
    writer[writer_num].Put<int>(init_int_const, vINT.data());
    writer[writer_num].Put<double>(init_double_const, vDOUBLE.data());
    writer[writer_num].EndStep();
    writer[writer_num].BeginStep();
    writer[writer_num].Put<double>(x, xml_in);
    writer[writer_num].Put<double>(y, yml_in);
    writer[writer_num].Put<double>(z, zml_in);
    writer[writer_num].Put<double>(p, pr);
    writer[writer_num].Put<double>(vx, v);
    writer[writer_num].Put<double>(vy, u);
    writer[writer_num].Put<double>(vz, w);
    writer[writer_num].Put<double>(t, tem);
    writer[writer_num].Put<double>(bm1, bm);
    writer[writer_num].Put<int>(InsituCounter, vInSituCounter.data());
    writer[writer_num].EndStep();
    /*
    adios2::Variable<double> x_xyz = io_xyz.DefineVariable<double>("X", {total1}, {start1}, {count1});
    adios2::Variable<double> y_xyz = io_xyz.DefineVariable<double>("Y", {total1}, {start1}, {count1});
    adios2::Variable<double> z_xyz = io_xyz.DefineVariable<double>("Z", {total1}, {start1}, {count1});
    writer0.Put<double>(x_xyz, xml_in);
    writer0.Put<double>(y_xyz, yml_in);
    writer0.Put<double>(z_xyz, zml_in);
    writer0.EndStep();
    writer0.Close();
    */
    firstStep=true;
    ++writer_num;
    ++insituCount;
    if(!rank) std::cout << "In-Situ setting done" << std::endl;
    std::cout << "Nek rank: " << rank << " count: " << nelt << " , start: " << nelbt << " , total: " << nelgt << " , nelbv: " << *nelbt_in << std::endl;
}


extern "C" void adios2_update_(
    const double *xml_in,
    const double *yml_in,
    const double *zml_in,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *temp,
    const double *bm
){  
    startT = std::clock();
    writer[insituCount%writer_num].BeginStep();
    vInSituCounter[0] = insituCount;

    /*
    if(firstStep){
    	writer.Put<double>(x, xml_in,adios2::Mode::Sync);
    	writer.Put<double>(y, yml_in,adios2::Mode::Sync);
    	writer.Put<double>(z, zml_in,adios2::Mode::Sync);
	firstStep=false;
    }
    */
    writer[insituCount%writer_num].Put<double>(p, pr);
    writer[insituCount%writer_num].Put<double>(vx, v);
    writer[insituCount%writer_num].Put<double>(vy, u);
    writer[insituCount%writer_num].Put<double>(vz, w);
    writer[insituCount%writer_num].Put<double>(t, temp);
    writer[insituCount%writer_num].Put<double>(bm1, bm);
    writer[insituCount%writer_num].Put<int>(InsituCounter, vInSituCounter.data());
    writer[insituCount%writer_num].EndStep();
    ++insituCount;
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_finalize_(){
    int closeIdx;
    for(closeIdx=0;closeIdx<writer_num;++closeIdx)
        writer[closeIdx].Close();
    
    std::cout <<  "rank: " << rank << " sin-situ time: " << dataTime << "s, total time: " << (std::clock() - startTotal) / (double) CLOCKS_PER_SEC << "s. " << std::endl;
}




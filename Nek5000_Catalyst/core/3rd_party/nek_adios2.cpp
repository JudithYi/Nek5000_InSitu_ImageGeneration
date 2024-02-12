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
std::vector<std::vector<adios2::Engine>> writer;
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
std::vector<int> writer_num;
std::vector<int> insituCount;
std::vector<int> stepCount;
std::vector<int> world_size;
std::vector<int> insitu_size;
std::vector<int> insitu_fre;
int adios_init;
int type_num = 0;

/* variables for reading */
adios2::IO ior;
adios2::Variable<double> vxr;
adios2::Variable<double> vyr;
adios2::Variable<double> vzr;
adios2::Variable<double> prr;
adios2::Variable<int> lglelr;
std::vector<double> vVXr;
std::vector<double> vVYr;
std::vector<double> vVZr;
std::vector<double> vVPrr;
std::vector<int> vVLGLELr;
int ifile = 0;

void init_multiple_type(const int type_num_in){
    type_num = type_num_in;
    writer_num.resize(type_num);
    insituCount.resize(type_num);
    world_size.resize(type_num);
    insitu_size.resize(type_num);
    insitu_fre.resize(type_num);
    stepCount.resize(type_num);
    writer.resize(type_num);
    for(i=0;i<type_num;++i){
        writer_num[i] = 0;
        insituCount[i] = 0;
        stepCount[i] = 0;
        insitu_fre[i] = 0;
    }
    adios_init = 0;
}

void adios_writer_init(
    MPI_Comm & comm,
    MPI_Comm & worldComm, 
    std::string engineName,
    const int iostep_in,
    const int type_index
){  
    if(adios_init == 0){
	startTotal = std::clock();
	std::string configFile="config/config.xml";
        adios = adios2::ADIOS(configFile, comm);
    	io = adios.DeclareIO("writer");
        MPI_Comm_rank(comm, &rank);
        MPI_Comm_size(comm, &size);
    }
    if(writer_num[type_index] == 0){
        MPI_Comm_size(worldComm, &world_size[type_index]);
        // int	worldComm_f = MPI_Comm_c2f(worldComm);
        // int comm_f = MPI_Comm_c2f(comm);
        // std::cout << "From writer[" << writer_num << "]: " << comm_f << " " << worldComm_f << " " << engineName << std::endl;
        insitu_size[type_index] = world_size[type_index] - size;
        insitu_fre[type_index] = iostep_in;
    }
    writer[type_index].resize(writer_num[type_index]+1);
    writer[type_index][writer_num[type_index]] = io.Open(engineName, adios2::Mode::Write,comm);
    //if(io.EngineType()=="InSituMPI")
    io.set_worldComm(&writer[type_index][writer_num[type_index]],&worldComm);
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
    const int *type_idx
){
    //writer0.BeginStep();
    int type_index = *type_idx;
    if(adios_init == 0){
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
        adios_init = 1;
    }
    vInSituCounter[0] = 0;

    adios2::StepStatus status = writer[type_index][writer_num[type_index]].BeginStep();
    if (status != adios2::StepStatus::OK){
	    std::cout << "Writer has not done begin step." << std::endl;
        writer[type_index][writer_num[type_index]].EndStep();
        writer[type_index][writer_num[type_index]].Close();
        //MPI_Finalize();
        return;
    }

    writer[type_index][writer_num[type_index]].Put<int>(init_int_const, vINT.data());
    writer[type_index][writer_num[type_index]].Put<double>(init_double_const, vDOUBLE.data());
    writer[type_index][writer_num[type_index]].EndStep();
    writer[type_index][writer_num[type_index]].BeginStep();
    writer[type_index][writer_num[type_index]].Put<double>(x, xml_in);
    writer[type_index][writer_num[type_index]].Put<double>(y, yml_in);
    writer[type_index][writer_num[type_index]].Put<double>(z, zml_in);
    writer[type_index][writer_num[type_index]].Put<double>(p, pr);
    writer[type_index][writer_num[type_index]].Put<double>(vx, v);
    writer[type_index][writer_num[type_index]].Put<double>(vy, u);
    writer[type_index][writer_num[type_index]].Put<double>(vz, w);
    writer[type_index][writer_num[type_index]].Put<double>(t, tem);
    writer[type_index][writer_num[type_index]].Put<double>(bm1, bm);
    writer[type_index][writer_num[type_index]].Put<int>(InsituCounter, vInSituCounter.data());
    writer[type_index][writer_num[type_index]].EndStep();

    ++writer_num[type_index];
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
    const double *bm,
    const int* type_idx
){  
    startT = std::clock();
    int type_index = *type_idx;
    ++stepCount[type_index];
    if(stepCount[type_index]%insitu_fre[type_index]!=0){
        dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
        return;
    }
    ++insituCount[type_index];
    int idx_writer = insituCount[type_index]%writer_num[type_index];
    writer[type_index][idx_writer].BeginStep();
    vInSituCounter[0] = stepCount[type_index];
    writer[type_index][idx_writer].Put<double>(p, pr);
    writer[type_index][idx_writer].Put<double>(vx, v);
    writer[type_index][idx_writer].Put<double>(vy, u);
    writer[type_index][idx_writer].Put<double>(vz, w);
    writer[type_index][idx_writer].Put<double>(t, temp);
    writer[type_index][idx_writer].Put<double>(bm1, bm);
    writer[type_index][idx_writer].Put<int>(InsituCounter, vInSituCounter.data());
    writer[type_index][idx_writer].EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

extern "C" void adios2_finalize_(){
    int closeIdx, typeIdx;
    for(typeIdx = 0; typeIdx < type_num; ++typeIdx){
        for(closeIdx = 0; closeIdx < writer_num[typeIdx]; ++closeIdx)
            writer[typeIdx][closeIdx].Close();
    }
    std::cout <<  "rank: " << rank << " sin-situ time: " << dataTime << "s, total time: " << (std::clock() - startTotal) / (double) CLOCKS_PER_SEC << "s. " << std::endl;
}

/* Here are the adios2 for lossy compression part */

extern "C" void adios2_update_lossy_(
    const int *lglel,
    const double *pr,
    const double *v,
    const double *u,
    const double *w,
    const double *temp,
    const int* type_idx
){
    startT = std::clock();
    int type_index = * type_idx;
    ++insituCount[type_index];
    int idx_writer = insituCount[type_index]%writer_num[type_index];
    writer[type_index][idx_writer].BeginStep();
    vInSituCounter[0] = insituCount[type_index] * insitu_fre[type_index];

    writer[type_index][idx_writer].Put<double>(p, pr);
    writer[type_index][idx_writer].Put<double>(vx, v);
    writer[type_index][idx_writer].Put<double>(vy, u);
    writer[type_index][idx_writer].Put<double>(vz, w);
    writer[type_index][idx_writer].Put<double>(t, temp);
    writer[type_index][idx_writer].Put<int>(InsituCounter, vInSituCounter.data());
    writer[type_index][idx_writer].EndStep();
    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
}

/* functions for reading */
extern "C" void adios2_read_(
    int *lglelrr,
    double *pr,
    double *v,
    double *u,
    double *w,
    //const double *temp,
    const int *nval,
    const int *nelvin,
    const int *nelb,
    const int *nelgv,
    const int *nelgt,
    const int *comm_int,
    char *fname
){
    startT = std::clock();
    std::string configFile="config/config.xml";
    MPI_Comm comm = MPI_Comm_f2c(*comm_int);
    //adios = adios2::ADIOS(configFile, comm);
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    // In config, writer is opening the insituMPI engine.
    // I guess this is due to the fact that this is just a write
    // to let the compressor executable access the data.
//    Coment all concerning ior since it was declared in the setup
//    ior = adios.DeclareIO("inputReader");
//    // In config, writer 0 uses the BPfile engine, wich i guess stands for//
//    // bzip2
//    if (!ior.InConfigFile())
//    {
//        // if not defined by user, we can change the default settings
//        // BPFile is the default engine
//        ior.SetEngine("BPFile");
//        ior.SetParameters({{"num_threads", "1"}});
//
//        // ISO-POSIX file output is the default transport (called "File")
//        // Passing parameters to the transport
//    }
    // Understand how the elements are divided among ranks.
    // set the pointers to the starting points in vectors of for each rank.
    //unsigned int nelv = static_cast<unsigned int>((*nelgv) / size);
    ior = adios.DeclareIO("inputReader");
    // In config, writer 0 uses the BPfile engine, wich i guess stands for
    // bzip2
    if (!ior.InConfigFile())
    {
        // if not defined by user, we can change the default settings
        // BPFile is the default engine
        ior.SetEngine("BPFile");
        ior.SetParameters({{"num_threads", "1"}});

        // ISO-POSIX file output is the default transport (called "File")
        // Passing parameters to the transport
    }
    unsigned int nelv = static_cast<unsigned int>((*nelvin));
    unsigned int start = static_cast<unsigned int>(*nelb);
    //if((*nelgv)%size != 0){
    //    if(rank < ((*nelgv) % size)){
    //        ++nelv;
    //        start += static_cast<unsigned int>(rank);
    //    }else{
    //        start += static_cast<unsigned int>((*nelgv)%size);
    //    }
    //}
    start *= static_cast<unsigned int>(*nval);
    // n is count, i.e number of entries in the array in my rank
    unsigned int n = static_cast<unsigned int> (*nval) * nelv;
    // gn is the total size of the arrays, not per io rank
    unsigned int gn = static_cast<unsigned int>((*nelgv)*(*nval));

    // These are the indicators for the global indices
    unsigned int nelv3 = static_cast<unsigned int>((*nelvin));
    unsigned int start3 = static_cast<unsigned int>(*nelb);
    unsigned int n3 = static_cast<unsigned int> (nelv);
    unsigned int gn3 = static_cast<unsigned int>((*nelgv));

    if (rank==0){

    std::cout << "what adios is getting:"<< std::endl;
    std::cout << "nvals:"<< static_cast<unsigned int> (*nval) << std::endl;
    std::cout << "nelgv:"<< static_cast<unsigned int> (*nelgv) << std::endl;
    std::cout << "nelgt:"<< static_cast<unsigned int> (*nelgt) << std::endl;

    std::cout << "what adios is calculating:"<< std::endl;
    std::cout << "nelv:"<< nelv << std::endl;
    std::cout << "start:"<< start << std::endl;
    std::cout << "count:"<< n << std::endl;
    std::cout << "total number of entries:"<< gn << std::endl;

    }
    //std::cout << "rank"<< rank << ": " << gn << ", " << start << "," << n << std::endl;

    //std::cout << "Reading the file " << std::endl;

    // This file is here just to correctly advance the reader
    // I have an error if I send the name from fortram. I would need to trim the string.
    ifile=ifile+1;
    // Make a new name every time the reader is called
    std:: string fileName= "out.f0000"+ std::to_string(ifile) +".bp";
    // Read the new name
    adios2::Engine readr = ior.Open(fileName,adios2::Mode::Read,comm);

    //adios2::Engine readr = ior.Open("out.f00001.bp",adios2::Mode::Read,comm);
    //adios2::Engine readr = ior.Open(fname,adios2::Mode::Read,comm);
    int step = 1;
    bool firstStep = true;

    // I am commenting the while loop and all the breaks since they
    // appear to create problems when reading multiple files in succession.
    // the problem is somewhere in the "steps"

    //while (true)
    //{
     	//adios2::StepStatus status =
        //readr.BeginStep(adios2::StepMode::Read);
    //    if (status != adios2::StepStatus::OK)
    //    {
    //        break;
    //    }

        vxr = ior.InquireVariable<double>("VX_OUT");
        if (!vxr)
        {
            std::cout << "Error: NO variable VX_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        vyr = ior.InquireVariable<double>("VY_OUT");
        if (!vyr)
        {
            std::cout << "Error: NO variable VY_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        vzr = ior.InquireVariable<double>("VZ_OUT");
        if (!vzr)
        {
            std::cout << "Error: NO variable VZ_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        prr = ior.InquireVariable<double>("P_OUT");
        if (!prr)
        {
            std::cout << "Error: NO variable P_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }
        lglelr = ior.InquireVariable<int>("LGLEL_OUT");
        if (!lglelr)
        {
            std::cout << "Error: NO variable LGLEL_OUT found. Unable to proceed. " << std::endl;
    //        break;
        }

        if (firstStep)
        {
        	// Promise that we are not going to change the variable sizes
                // nor add new variables
                readr.LockReaderSelections();

		vVXr.resize(n);
		vVYr.resize(n);
		vVZr.resize(n);
		vVPrr.resize(n);
		vVLGLELr.resize(n3);

		firstStep=false;
        }


        // std::cout << "Select how much is yours "
           //               << std::endl;
        vxr.SetSelection({{start}, {n}});
        vyr.SetSelection({{start}, {n}});
        vzr.SetSelection({{start}, {n}});
        prr.SetSelection({{start}, {n}});
        lglelr.SetSelection({{start3}, {n3}});


        // std::cout << "get the data "
           //              << std::endl;
        readr.Get<double>(vxr,vVXr.data());
        readr.Get<double>(vyr,vVYr.data());
        readr.Get<double>(vzr,vVZr.data());
        readr.Get<double>(prr,vVPrr.data());
        readr.Get<int>(lglelr,vVLGLELr.data());
        //readr.Get<double>(vzr,w);

        //readr.EndStep();
        //++step;
    //}
    readr.Close();
    //++step;

    //copy the data into w
    if (rank==0){
    std::cout <<"Copying data into nek vector " << std::endl;
    }
    //if (rank==0){
    for (int i=0; i<vVZr.size(); i++){
    //    std::cout <<"entry "<< i <<"is " << vVZr[i] << " ";
    //    std::cout <<"size of read vector " <<"is " << vvzr.size() << std::endl;
    pr[i]=vVPrr[i];
    v[i]=vVXr[i];
    u[i]=vVYr[i];
    w[i]=vVZr[i];
    }
    for (int i=0; i<vVLGLELr.size(); i++){
    //    std::cout <<"entry "<< i <<"is " << vVZr[i] << " ";
    //    std::cout <<"size of read vector " <<"is " << vvzr.size() << std::endl;
    lglelrr[i]=vVLGLELr[i];
    }
    //}

    dataTime += (std::clock() - startT) / (double) CLOCKS_PER_SEC;
    //}

}

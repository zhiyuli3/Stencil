# Stencil

## Some Useful Scripts
module load  
python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm
## MPI Function Usage
`int MPI_Init(int *argc, char*argv[]);`  
这句应该放在所有并行操作的最前面。  
`int MPI_Initialized(int *flag);`  
是否初始化， 结果保存在flag中  
`MPI_Comm_rank( MPI_Comm communicator,int* rank)`
这个函数会返回 communicator 中当前进程的 rank。 communicator 中每个进程会以此得到一个从0开始递增的数字作为 rank 值。rank 值主要是用来指定发送或者接受信息时对应的进程。size是总共到进程数目。
`MPI_Comm_size`
总的进程数目
`MPI_Send(void* data,int count,MPI_Datatype datatype,int destination,int tag,MPI_Comm communicator)`
Example:MPI_Send(message,strlen(message)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD);
`MPI_Recv(void* data,int count,MPI_Datatype datatype,int source,int tag,MPI_Comm communicator,MPI_Status* status)`
`int MPI_Abort(MPI_Comm comm, int errorcode)`
终止MPI环境及MPI程序的执行
`MPI_Wtime()`
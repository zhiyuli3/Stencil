# Stencil

## Some Useful Scripts
module load  
python check.py --ref-stencil-file stencil_1024_1024_100.pgm --stencil-file stencil.pgm
## MPI Function Usage
`int MPI_Init(int *argc, char*argv[]);`  
这句应该放在所有并行操作的最前面。 
_______
`int MPI_Initialized(int *flag);`  
是否初始化， 结果保存在flag中 
_______ 
`MPI_Comm_rank( MPI_Comm communicator,int* rank)`    
这个函数会返回 communicator 中当前进程的 rank。 communicator 中每个进程会以此得到一个从0开始递增的数字作为 rank 值。rank 值主要是用来指定发送或者接受信息时对应的进程。size是总共到进程数目。  
_______
`MPI_Comm_size`  
总的进程数目  
_______
`MPI_Send(void* data,int count,MPI_Datatype datatype,int destination,int tag,MPI_Comm communicator)`  
Example:MPI_Send(message,strlen(message)+1, MPI_CHAR, dest, tag, MPI_COMM_WORLD); 
_______ 
`MPI_Recv(void* data,int count,MPI_Datatype datatype,int source,int tag,MPI_Comm communicator,MPI_Status* status)`  
_______
`int MPI_Abort(MPI_Comm comm, int errorcode)`  
终止MPI环境及MPI程序的执行  
_______
`MPI_Wtime()` 
_______ 
`int MPI_Isend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag,MPI_Comm comm, MPI_Request *request)`  
输入参数：  
buf：发送缓冲区的首地址   
count：需要发送的字节数  
datatype：每个发送元素的数据类型  
dest：目标的rank（id）   
tag：消息标识（integer）  
comm：通信域  
输出参数：  
request:communication request (handle)  
_______
`int MPI_Irecv(void *buf, int count, MPI_Datatype datatype, int source,int tag, MPI_Comm comm, MPI_Request *request)`  
输出参数：  
buf：接收缓冲区的首地址  
count：接收缓冲区存放字节数(integer)  
datatype：每个接收元素的数据类型  
source：发送者的rank (integer)  
tag：消息标识（integer）  
comm：通信域  
输出参数：  
request：communication request (handle)  
`MPI_Wait`  
_______
`int MPI_Wait(MPI_Request *request, MPI_Status *status)`   
_______
`int MPI_Test(MPI_Request *request, int *flag, MPI_Status *status)`   
_______
`int MPI_Ssend(const void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)`  
Blocking synchronous send  
当一个匹配到一个接收者并且开始接收时，同步模式就算完成，但是他的完成并不代表接收者已经完成接收，可以释放缓存，而且如果接收和发送都是堵塞的，这个可以提供同步语义
________  
```
MPI_Sendrecv( void *sendbuf //initial address of send buffer   
int sendcount //number of entries to send  
MPI_Datatype sendtype //type of entries in send buffer  
int dest //rank of destination
int sendtag //send tag
void *recvbuf //initial address of receive buffer
int recvcount //max number of entries to receive  
MPI_Datatype recvtype //type of entries in receive buffer 
```  
_________   
`int MPI_Waitall(int count, MPI_Request array_of_requests[], MPI_Status array_of_statuses[])`   
Waits for all given MPI Requests to complete   
## Input Parameters:  
count list length (integer)  
array_of_requests array of request handles (array of handles)  

## Output Parameters:  
array_of_statuses
array of status objects (array of Statuses). May be MPI_STATUSES_IGNORE.  
__________
`MPI_Win_create(*base, total_bytes, unit_bytes, info, comm, *win_handle)`  
创建一个 Window 并暴露一段连续的内存给 MPI   
单边MPI传递
https://enigmahuang.github.io/2017/06/26/MPI3-OSC/   

`MPI_Bcast(void *buffer, int count, MPI_Datatype datatype, int root, MPI_Comm comm)`
Broadcasts a message from the process with rank "root" to all other processes of the communicator
________________________________________
`MPI_Wtick()`  
Returns the resolution of MPI_Wtime
_______________________________________    
`MPI_Reduce(const void *sendbuf, void *recvbuf, int count, MPI_Datatype datatype,MPI_Op op, int root, MPI_Comm comm)`  
______________________________________  
```
MPI_SCATTER(sendbuf,sendcount,sendtype,recvbuf,recvcount,recvtype,
            root,comm)
 IN   sendbuf     发送消息缓冲区的起始地址(可变,仅对于根进程)
 IN   sendcount   发送到各个进程的数据个数(整型,仅对于根进程)
 IN   sendtype    发送消息缓冲区中的数据类型(句柄,仅对于根进程)
 OUT  recvbuf     接收消息缓冲区的起始地址(可变)
 IN   recvcount   待接收的元素个数(整型)
 IN   recvtype    接收元素的数据类型(句柄)
 IN   root        发送进程的序列号(整型)
 IN   comm        通信子(句柄)
int MPI_Scatter(void* sendbuf, int sendcount, MPI_Datatype sendtype,
                void* recvbuf, int recvcount, MPI_Datatype recvtype,
                int root, MPI_Comm comm)  
```  
_________________________  
```  
MPI_GATHER(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype, root , comm)
　IN　sendbuf   　发送消息缓冲区的起始地址(可变)
　IN　sendcount 　发送消息缓冲区中的数据个数(整型)
　IN　sendtype　  发送消息缓冲区中的数据类型(句柄) 
　OUT recvbuf 　　接收消息缓冲区的起始地址(可变,仅对于根进程) 
　IN　recvcount 　待接收的元素个数(整型,仅对于根进程)
　IN　recvtype 　 接收元素的数据类型(句柄,仅对于根进程)
　IN　root　　　   接收进程的序列号(整型)
　IN　comm 　　 　 通信子(句柄)
int MPI_Gather(void* sendbuf, int sendcount, MPI_Datatype sendtype, 
               void* recvbuf, int recvcount, MPI_Datatype recvtype, 
               int root, MPI_Comm comm)
```  
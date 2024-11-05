
#include "util.h"



struct ringbuffer{
	char* buffer;	//缓冲内存
	int size;		//缓冲区大小
	int front;		//读者索引，即消费者索引
	int rear;		//写者索引，即生产者线程
	pthread_mutex_t  mt;  //互斥量
	pthread_cond_t  cond;	//条件变量
};




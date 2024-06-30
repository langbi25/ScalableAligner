#include "PthreadPool.h"

PthreadPool::PthreadPool()
{
    thread_num = 0;
    running_num = 0;
    shutdown = 0;
}

PthreadPool::~PthreadPool()
{
    pthread_mutex_lock(&lock);    // 先上锁， 防止有程序占用
    pthread_mutex_destroy(&lock); // 销毁
    pthread_cond_destroy(&notify);
    delete[] threads;
}

int PthreadPool::Init(unsigned int num)
{
    // 初始化互斥锁和条件变量
    do
    {
        if (num <= 0)
            break;
        if (pthread_mutex_init(&lock, NULL))
            break;
        if (pthread_cond_init(&notify, NULL))
            break;

        // 初始化线程数组
        threads = new pthread_t[num];

        // 创建线程
        for (int i = 0; i < num; i++)
        {
            if (pthread_create(threads + i, NULL, threadpool_thread, (void *)this) != 0)
            {
                // 创建不成功则销毁
                Destory();
                break;
            }
            running_num++;
            thread_num++;
        }
        return 0; // 成功

    } while (0);
    thread_num = 0;
    return Pthreadpool_invalid;
}

int PthreadPool::Destory(PthreadPool_Shutdown flag)
{
    do
    {
        // 取得互斥锁资源
        if (pthread_mutex_lock(&lock) != 0)
        {
            return Pthreadpool_lock_failure;
        }

        shutdown = flag; // 标记标记

        /* 唤醒所有因条件变量阻塞的线程，并释放互斥锁 */
        if ((pthread_cond_broadcast(&notify) != 0) || (pthread_mutex_unlock(&lock) != 0))
            break;

        /* 等待所有线程结束 */
        for (int i = 0; i < thread_num; i++)
            if (pthread_join(threads[i], NULL) != 0)
                break;
        return 0;
    } while (0);
    return -1;
}

int PthreadPool::AddTask(void (*function)(void *), void *argument)
{

    if (thread_num == 0 || function == NULL)
        return Pthreadpool_invalid;
    /* 必须先取得互斥锁所有权 */
    if (pthread_mutex_lock(&lock) != 0)
        return Pthreadpool_lock_failure;

    // 检查是否关闭了线程池
    if (shutdown){
        return Pthreadpool_shutdown;
    }

    // 新加入
    Pthreadpool_Runable newRunable;
    newRunable.function = function;
    newRunable.argument = argument;
    // 加入队列
    process_num++;
    thread_queue.push(newRunable);
    // 发出signal
    if (pthread_cond_signal(&notify) != 0)
        return Pthreadpool_lock_failure;
    pthread_mutex_unlock(&lock);
    return 0;
}

// 线程运行函数
void *PthreadPool::threadpool_thread(void *threadpool)
{
    PthreadPool *pool = (PthreadPool *)threadpool; // 获取当前实例
    while (1){
        /* 取得互斥锁资源 */
        pthread_mutex_lock(&(pool->lock));
        while ((pool->thread_queue.empty()) && (!pool->shutdown)){
            /* 任务队列为空，且线程池没有关闭时阻塞在这里 */
            pthread_cond_wait(&(pool->notify), &(pool->lock));
        }

        /* 关闭的处理 */
        if((pool->shutdown == immediate_shutdown) ||
           ((pool->shutdown == graceful_shutdown) && (pool->thread_queue.empty()))) {
            break;
        }

        // 取队列中的任务
        Pthreadpool_Runable Runable;
        if (!pool->thread_queue.empty()){
            Runable = pool->thread_queue.front();
                        // fprintf(stderr, "pool.thread_queue.empty():\n");

            pool->thread_queue.pop(); // 出队
        }

        // running_num
        pool->waiting_num--;
        /* 释放互斥锁 */
        pthread_mutex_unlock(&(pool->lock));

        // 开始运行任务
        (*(Runable.function))(Runable.argument);
        // pool->process_num--;
        pool->process_num--;
        // pool->process_num--;


        // 结束，回到等待
    }
    // 更新正在运行的线程数
    pool->running_num--;

    /* 释放互斥锁 */
    pthread_mutex_unlock(&(pool->lock));
    pthread_exit(NULL);
    return (NULL);
}
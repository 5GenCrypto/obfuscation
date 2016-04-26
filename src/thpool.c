/* ********************************
 * Author:       Johan Hanssen Seferidis
 * License:	     MIT
 * Description:  Library providing a threading pool where you can add
 *               work. For usage, check the thpool.h file or README.md
 *
 *//** @file thpool.h *//*
 * 
 ********************************/

#include <unistd.h>
#include <signal.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <errno.h>
#include <time.h> 
#if defined(__linux__)
#include <sys/prctl.h>
#endif
#include "thpool.h"

#ifdef THPOOL_DEBUG
#define THPOOL_DEBUG 1
#else
#define THPOOL_DEBUG 0
#endif

static volatile int threads_keepalive;
static volatile int threads_on_hold;


/* ========================== STRUCTURES ============================ */


/* Binary semaphore */
typedef struct bsem {
	pthread_mutex_t mutex;
	pthread_cond_t   cond;
	int v;
} bsem;


/* Job */
typedef struct job {
	struct job*  prev;                   /* pointer to previous job   */
	void*  (*function)(void* arg);       /* function pointer          */
	void*  arg;                          /* function's argument       */
    char *tag;
} job;


/* Job queue */
typedef struct jobqueue{
	pthread_mutex_t rwmutex;             /* used for queue r/w access */
	job  *front;                         /* pointer to front of queue */
	job  *rear;                          /* pointer to rear  of queue */
	bsem *has_jobs;                      /* flag as binary semaphore  */
	int   len;                           /* number of jobs in queue   */
} jobqueue;


/* Thread */
typedef struct thread{
	int       id;                        /* friendly id               */
	pthread_t pthread;                   /* pointer to actual thread  */
	struct thpool_* thpool_p;            /* access to thpool          */
} thread;

typedef struct tag {
    char *name;
    int len;
    void * (*function)(void *arg);
    void *arg;
} tag_t;

typedef struct taglist_ {
    int num;
    pthread_mutex_t lock;
    tag_t *tags;
} taglist_t;

// XXX: We don't clear the tag list, so if we have more than 1024 tags all hell
// breaks loose!
#define TAGLIST_SIZE 1024

/* Threadpool */
typedef struct thpool_{
	thread**   threads;                  /* pointer to threads        */
	volatile int num_threads_alive;      /* threads currently alive   */
	volatile int num_threads_working;    /* threads currently working */
	pthread_mutex_t  thcount_lock;       /* used for thread count etc */
	pthread_cond_t  threads_all_idle;    /* signal to thpool_wait     */
	jobqueue*  jobqueue_p;               /* pointer to the job queue  */
    taglist_t* tlist;
} thpool_;


/* ========================== PROTOTYPES ============================ */


static int  thread_init(thpool_* thpool_p, struct thread** thread_p, int id);
static void* thread_do(struct thread* thread_p);
static void  thread_hold();
static void  thread_destroy(struct thread* thread_p);

static int   jobqueue_init(thpool_* thpool_p);
static void  jobqueue_clear(thpool_* thpool_p);
static void  jobqueue_push(thpool_* thpool_p, struct job* newjob_p);
static struct job* jobqueue_pull(thpool_* thpool_p);
static void  jobqueue_destroy(thpool_* thpool_p);

static void  bsem_init(struct bsem *bsem_p, int value);
static void  bsem_reset(struct bsem *bsem_p);
static void  bsem_post(struct bsem *bsem_p);
static void  bsem_post_all(struct bsem *bsem_p);
static void  bsem_wait(struct bsem *bsem_p);


/* ========================== THREADPOOL ============================ */


/* Initialise thread pool */
struct thpool_*
thpool_init(int num_threads)
{
	threads_on_hold   = 0;
	threads_keepalive = 1;

	if (num_threads < 0){
		num_threads = 0;
	}

	/* Make new thread pool */
	thpool_* thpool_p;
	thpool_p = (struct thpool_*)malloc(sizeof(struct thpool_));
	if (thpool_p == NULL) {
		fprintf(stderr, "thpool_init(): Could not allocate memory for thread pool\n");
		return NULL;
	}
	thpool_p->num_threads_alive   = 0;
	thpool_p->num_threads_working = 0;

	/* Initialise the job queue */
	if (jobqueue_init(thpool_p) == -1) {
		fprintf(stderr, "thpool_init(): Could not allocate memory for job queue\n");
		free(thpool_p);
		return NULL;
	}

	/* Make threads in pool */
	thpool_p->threads = (struct thread**)malloc(num_threads * sizeof(struct thread *));
	if (thpool_p->threads == NULL){
		fprintf(stderr, "thpool_init(): Could not allocate memory for threads\n");
		jobqueue_destroy(thpool_p);
		free(thpool_p->jobqueue_p);
		free(thpool_p);
		return NULL;
	}

    /* Make tag list */
    taglist_t *tlist;
    tlist = (taglist_t *) malloc(sizeof(taglist_t));
    if (tlist == NULL) {
        fprintf(stderr, "thpool_init(): Could not allocate memory for job queue\n");
		jobqueue_destroy(thpool_p);
		free(thpool_p->jobqueue_p);
        free(thpool_p->threads);
		free(thpool_p);
        return NULL;
    }
    tlist->num = TAGLIST_SIZE;
    tlist->tags = (tag_t *) calloc(tlist->num, sizeof(tag_t));
    if (tlist->tags == NULL) {
        fprintf(stderr, "thpool_init(): Could not allocate memory for tag list\n");
		jobqueue_destroy(thpool_p);
		free(thpool_p->jobqueue_p);
        free(thpool_p->threads);
		free(thpool_p);
        free(tlist);
        return NULL;
    }
    pthread_mutex_init(&tlist->lock, NULL);
    thpool_p->tlist = tlist;

	pthread_mutex_init(&(thpool_p->thcount_lock), NULL);
	pthread_cond_init(&thpool_p->threads_all_idle, NULL);
	
	/* Thread init */
	for (int n = 0; n < num_threads; n++) {
		thread_init(thpool_p, &thpool_p->threads[n], n);
		if (THPOOL_DEBUG)
			printf("THPOOL_DEBUG: Created thread %d in pool \n", n);
	}
	
	/* Wait for threads to initialize */
	while (thpool_p->num_threads_alive != num_threads) {}

	return thpool_p;
}

int
thpool_add_tag(thpool_ *thpool_p, char *tag, int length,
               void * (fn)(void *arg), void *arg)
{
    for (int i = 0; i < thpool_p->tlist->num; ++i) {
        struct tag *t = &thpool_p->tlist->tags[i];
        if (t->name == NULL) {
            t->name = (char *) calloc(strlen(tag) + 1, sizeof(char));
            (void) strcpy(t->name, tag);
            t->len = length;
            t->function = fn;
            t->arg = arg;
            return 0;
        } else if (strcmp(tag, t->name) == 0) {
            fprintf(stderr, "thpool_add_tag(): tag '%s' already in tag list\n", tag);
            return -1;
        }
    }

    fprintf(stderr, "thpool_add_tag(): tag list full!\n");
    return -1;
}


static int
tag_decrement(thpool_ *thpool_p, char *tag)
{
    for (int i = 0; i < thpool_p->tlist->num; ++i) {
        struct tag *t = &thpool_p->tlist->tags[i];
        if (strcmp(t->name, tag) == 0) {
            pthread_mutex_lock(&thpool_p->tlist->lock);
            t->len--;
            pthread_mutex_unlock(&thpool_p->tlist->lock);
            if (t->len == 0) {
                t->function(t->arg);
            }
            return 0;
        }
    }
    fprintf(stderr, "tag_decrement(): tag does not exist in tag list\n");
    return -1;
}


/* Add work to the thread pool */
int
thpool_add_work(thpool_* thpool_p, void *(*function_p)(void*), void* arg_p,
                char *tag)
{
	job *newjob;

	newjob = (struct job*) malloc(sizeof(struct job));
	if (newjob == NULL) {
		fprintf(stderr, "thpool_add_work(): Could not allocate memory for new job\n");
		return -1;
	}

	/* add function and argument */
	newjob->function = function_p;
	newjob->arg = arg_p;
    newjob->tag = (char *) calloc(strlen(tag) + 1, sizeof(char));
    (void) strcpy(newjob->tag, tag);

	/* add job to queue */
	pthread_mutex_lock(&thpool_p->jobqueue_p->rwmutex);
	jobqueue_push(thpool_p, newjob);
	pthread_mutex_unlock(&thpool_p->jobqueue_p->rwmutex);

	return 0;
}


/* Wait until all jobs have finished */
void thpool_wait(thpool_* thpool_p){
	pthread_mutex_lock(&thpool_p->thcount_lock);
	while (thpool_p->jobqueue_p->len || thpool_p->num_threads_working) {
		pthread_cond_wait(&thpool_p->threads_all_idle, &thpool_p->thcount_lock);
	}
	pthread_mutex_unlock(&thpool_p->thcount_lock);
}


/* Destroy the threadpool */
void thpool_destroy(thpool_* thpool_p){
	/* No need to destory if it's NULL */
	if (thpool_p == NULL) return ;

	volatile int threads_total = thpool_p->num_threads_alive;

	/* End each thread's infinite loop */
	threads_keepalive = 0;
	
	/* Give one second to kill idle threads */
	double TIMEOUT = 1.0;
	time_t start, end;
	double tpassed = 0.0;
	time (&start);
	while (tpassed < TIMEOUT && thpool_p->num_threads_alive){
		bsem_post_all(thpool_p->jobqueue_p->has_jobs);
		time (&end);
		tpassed = difftime(end,start);
	}
	
	/* Poll remaining threads */
	while (thpool_p->num_threads_alive){
		bsem_post_all(thpool_p->jobqueue_p->has_jobs);
		sleep(1);
	}

	/* Job queue cleanup */
	jobqueue_destroy(thpool_p);
	free(thpool_p->jobqueue_p);
	
	/* Deallocs */
	int n;
	for (n=0; n < threads_total; n++){
		thread_destroy(thpool_p->threads[n]);
	}
	free(thpool_p->threads);
	free(thpool_p);
}


/* Pause all threads in threadpool */
void
thpool_pause(thpool_* thpool_p)
{
	for (int n = 0; n < thpool_p->num_threads_alive; n++) {
		pthread_kill(thpool_p->threads[n]->pthread, SIGUSR1);
	}
}

/* Resume all threads in threadpool */
void
thpool_resume(thpool_* thpool_p)
{
	threads_on_hold = 0;
}


/* ============================ THREAD ============================== */


/* Initialize a thread in the thread pool
 * 
 * @param thread        address to the pointer of the thread to be created
 * @param id            id to be given to the thread
 * @return 0 on success, -1 otherwise.
 */
static int thread_init (thpool_* thpool_p, struct thread** thread_p, int id){
	
	*thread_p = (struct thread*)malloc(sizeof(struct thread));
	if (thread_p == NULL){
		fprintf(stderr, "thpool_init(): Could not allocate memory for thread\n");
		return -1;
	}

	(*thread_p)->thpool_p = thpool_p;
	(*thread_p)->id       = id;

	pthread_create(&(*thread_p)->pthread, NULL, thread_do, (*thread_p));
	pthread_detach((*thread_p)->pthread);
	return 0;
}


/* Sets the calling thread on hold */
static void
thread_hold (int num)
{
	threads_on_hold = 1;
	while (threads_on_hold) {
		sleep(1);
	}
}


/* What each thread is doing
* 
* In principle this is an endless loop. The only time this loop gets interrupted
* is once thpool_destroy() is invoked or the program exits.
* 
* @param  thread        thread that will run this function
* @return nothing
*/
static void *
thread_do(struct thread* thread_p)
{
	/* Set thread name for profiling and debuging */
	char thread_name[128] = {0};
	sprintf(thread_name, "thread-pool-%d", thread_p->id);

#if defined(__linux__)
	/* Use prctl instead to prevent using _GNU_SOURCE flag and implicit declaration */
	prctl(PR_SET_NAME, thread_name);
#elif defined(__APPLE__) && defined(__MACH__)
	pthread_setname_np(thread_name);
#else
	fprintf(stderr, "thread_do(): pthread_setname_np is not supported on this system");
#endif

	/* Assure all threads have been created before starting serving */
	thpool_ *thpool_p = thread_p->thpool_p;
	
	/* Register signal handler */
	struct sigaction act;
	sigemptyset(&act.sa_mask);
	act.sa_flags = 0;
	act.sa_handler = thread_hold;
	if (sigaction(SIGUSR1, &act, NULL) == -1) {
		fprintf(stderr, "thread_do(): cannot handle SIGUSR1");
	}
	
	/* Mark thread as alive (initialized) */
	pthread_mutex_lock(&thpool_p->thcount_lock);
	thpool_p->num_threads_alive += 1;
	pthread_mutex_unlock(&thpool_p->thcount_lock);

	while (threads_keepalive) {

		bsem_wait(thpool_p->jobqueue_p->has_jobs);

		if (threads_keepalive) {

			job* job_p;
			
			pthread_mutex_lock(&thpool_p->thcount_lock);
			thpool_p->num_threads_working++;
			pthread_mutex_unlock(&thpool_p->thcount_lock);
			
			/* Read job from queue and execute it */
			pthread_mutex_lock(&thpool_p->jobqueue_p->rwmutex);
			job_p = jobqueue_pull(thpool_p);
			pthread_mutex_unlock(&thpool_p->jobqueue_p->rwmutex);
			if (job_p) {
				job_p->function(job_p->arg);
                tag_decrement(thpool_p, job_p->tag);
                free(job_p->tag);
				free(job_p);
			}
			
			pthread_mutex_lock(&thpool_p->thcount_lock);
			thpool_p->num_threads_working--;
			if (!thpool_p->num_threads_working) {
				pthread_cond_signal(&thpool_p->threads_all_idle);
			}
			pthread_mutex_unlock(&thpool_p->thcount_lock);

		}
	}
	pthread_mutex_lock(&thpool_p->thcount_lock);
	thpool_p->num_threads_alive --;
	pthread_mutex_unlock(&thpool_p->thcount_lock);

	return NULL;
}


/* Frees a thread  */
static void
thread_destroy (thread* thread_p)
{
	free(thread_p);
}


/* ============================ JOB QUEUE =========================== */


/* Initialize queue */
static int jobqueue_init(thpool_* thpool_p){
	
	thpool_p->jobqueue_p = (struct jobqueue*)malloc(sizeof(struct jobqueue));
	if (thpool_p->jobqueue_p == NULL){
		return -1;
	}
	thpool_p->jobqueue_p->len = 0;
	thpool_p->jobqueue_p->front = NULL;
	thpool_p->jobqueue_p->rear  = NULL;

	thpool_p->jobqueue_p->has_jobs = (struct bsem*)malloc(sizeof(struct bsem));
	if (thpool_p->jobqueue_p->has_jobs == NULL){
		return -1;
	}

	pthread_mutex_init(&(thpool_p->jobqueue_p->rwmutex), NULL);
	bsem_init(thpool_p->jobqueue_p->has_jobs, 0);

	return 0;
}


/* Clear the queue */
static void jobqueue_clear(thpool_* thpool_p){

	while(thpool_p->jobqueue_p->len){
		free(jobqueue_pull(thpool_p));
	}

	thpool_p->jobqueue_p->front = NULL;
	thpool_p->jobqueue_p->rear  = NULL;
	bsem_reset(thpool_p->jobqueue_p->has_jobs);
	thpool_p->jobqueue_p->len = 0;

}


/* Add (allocated) job to queue
 *
 * Notice: Caller MUST hold a mutex
 */
static void jobqueue_push(thpool_* thpool_p, struct job* newjob){

	newjob->prev = NULL;

	switch(thpool_p->jobqueue_p->len){

		case 0:  /* if no jobs in queue */
					thpool_p->jobqueue_p->front = newjob;
					thpool_p->jobqueue_p->rear  = newjob;
					break;

		default: /* if jobs in queue */
					thpool_p->jobqueue_p->rear->prev = newjob;
					thpool_p->jobqueue_p->rear = newjob;
					
	}
	thpool_p->jobqueue_p->len++;
	
	bsem_post(thpool_p->jobqueue_p->has_jobs);
}


/* Get first job from queue(removes it from queue)
 * 
 * Notice: Caller MUST hold a mutex
 */
static struct job *
jobqueue_pull(thpool_ *thpool_p)
{
	job *job_p;
	job_p = thpool_p->jobqueue_p->front;

	switch(thpool_p->jobqueue_p->len) {
		
    case 0:  /* if no jobs in queue */
        break;
		
    case 1:  /* if one job in queue */
        thpool_p->jobqueue_p->front = NULL;
        thpool_p->jobqueue_p->rear  = NULL;
        thpool_p->jobqueue_p->len = 0;
        break;
		
    default: /* if >1 jobs in queue */
        thpool_p->jobqueue_p->front = job_p->prev;
        thpool_p->jobqueue_p->len--;
        /* more than one job in queue -> post it */
        bsem_post(thpool_p->jobqueue_p->has_jobs);
					
	}
	
	return job_p;
}


/* Free all queue resources back to the system */
static void jobqueue_destroy(thpool_* thpool_p){
	jobqueue_clear(thpool_p);
	free(thpool_p->jobqueue_p->has_jobs);
}


/* ======================== SYNCHRONISATION ========================= */


/* Init semaphore to 1 or 0 */
static void bsem_init(bsem *bsem_p, int value) {
	if (value < 0 || value > 1) {
		fprintf(stderr, "bsem_init(): Binary semaphore can take only values 1 or 0");
		exit(1);
	}
	pthread_mutex_init(&(bsem_p->mutex), NULL);
	pthread_cond_init(&(bsem_p->cond), NULL);
	bsem_p->v = value;
}


/* Reset semaphore to 0 */
static void bsem_reset(bsem *bsem_p) {
	bsem_init(bsem_p, 0);
}


/* Post to at least one thread */
static void bsem_post(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_signal(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Post to all threads */
static void bsem_post_all(bsem *bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	bsem_p->v = 1;
	pthread_cond_broadcast(&bsem_p->cond);
	pthread_mutex_unlock(&bsem_p->mutex);
}


/* Wait on semaphore until semaphore has value 0 */
static void bsem_wait(bsem* bsem_p) {
	pthread_mutex_lock(&bsem_p->mutex);
	while (bsem_p->v != 1) {
		pthread_cond_wait(&bsem_p->cond, &bsem_p->mutex);
	}
	bsem_p->v = 0;
	pthread_mutex_unlock(&bsem_p->mutex);
}

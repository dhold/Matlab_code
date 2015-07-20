#ifndef Barrier_h
#define Barrier_h
#include <pthread.h>

typedef struct {
    int needed;
    int called;
    pthread_mutex_t mutex;
    pthread_cond_t cond;
} barrier_t;




int barrier_destroy(barrier_t *);
int barrier_init(barrier_t *, int);
int barrier_wait(barrier_t *);

#endif


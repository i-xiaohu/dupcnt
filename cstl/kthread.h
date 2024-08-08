//
// Created by ixiaohu on 2021/9/22.
//

#ifndef ZIP_SEEDING_KTHREAD_H
#define ZIP_SEEDING_KTHREAD_H

#ifdef __cplusplus
extern "C" {
#endif

void kt_for(int n_threads, void (*func)(void *, long, int), void *data, long n);

void kt_pipeline(int n_threads, void *(*func)(void *, int, void *), void *shared_data, int n_steps);

#ifdef __cplusplus
}
#endif

#endif //ZIP_SEEDING_KTHREAD_H

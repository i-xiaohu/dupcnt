//
// Created by ixiaohu on 2021/9/22.
//

#ifndef KOPEN_H
#define KOPEN_H

#ifdef __cplusplus
extern "C" {
#endif

	void *kopen(const char *fn, int *_fd);
	int kclose(void *a);

#ifdef __cplusplus
}
#endif

#endif //KOPEN_H

/*
 * list.h
 *
 *  Created on: Sep 6, 2016
 *      Author: Said Kerrache
 */

#ifndef HRGDS_H_
#define HRGDS_H_

#include <string>

class list {
public:
	int x;				// stored elementd in linked-list
	list* next;			// pointer to next elementd
	list();
	~list();
};

struct block {
	double x;
	int y;
};

struct ipair {
	int x;
	int y;
	short int t;
	std::string sp;
};

#endif

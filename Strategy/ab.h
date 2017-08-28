#ifndef AB_H_
#define AB_H_

void getPoint(const int , const int , const int* , int** , 
					 const int , const int , int, int, int& , int& );
int choose();
double ucth();
int get_first_expand(int);
void initial();
void iterate_imitate();
void update();
void expand(int node_pos);
int imitate(int player);
bool cmp(const int& A, const int& B);
void initial_keng();


#endif
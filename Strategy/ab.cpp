#include "Point.h"
#include "ab.h"
#include "judge.h"
#include <cstdlib>
#include <cstdio>
#include <ctime>
#include <algorithm>
#include <iostream>
#include <vector>
#include <fstream>
#ifndef PRINT
//#define PRINT	
#endif

#ifdef PRINT
#include <fstream>
#include <conio.h>
#include <atlstr.h>
#endif
using namespace std;
#ifdef PRINT	
ofstream file("a.in", ios::app);
ofstream kengdieba("b.in");
#endif

struct Node
{
	int player, round_sum, round_win, col;
	double tcb;
	//	��ǰ���֣���ǰ�ڵ���ʴ�������ǰ�ڵ��Լ�ʤ������, ����λ��
	std::vector<int> son;
	//	����
	
};

const double CONST_C = 0.9;
//double VAR_C;
//	����C: ����ϵ��	 EPS: utb��Ⱦ���
const int TREE_SIZE = 8001000;
//	����С	��������
const int ARRSIZE = 22;
//	�����С
Node tree[TREE_SIZE];
//	��
int size_pos, start_time, stack_top, rwin;
//	��λ��	��ʼʱ��	ջλ��	�ܳ���	ʤ��
int dad_stack[100000], next_expand[TREE_SIZE];
//	����ջ	����չ����һ��
int max_No, min_No, M, N, noX, noY, ans_x, ans_y;
//	max ��� min ��� M N noX noY �������
int top[ARRSIZE], Top[ARRSIZE], Board[ARRSIZE][ARRSIZE];
//	top;
int** board = NULL;
//	board;
int root, havechess = 10000, last_x, last_y;
//	��	������	������������λ��
int rand_seq[ARRSIZE], derta_seq[ARRSIZE];
int random_sequence[222], seq_point, seq_pos;
//	�������
//int keng_flag = 0;
//double MAXVAR = 2;




void initial_keng()
{
	//keng_flag = 0;
	//MAXVAR = 2;
	max_No = 0;
	min_No = 1;
	M = N = noX = noY = -1;
	int flag = 0;
	if (board == NULL)
		flag = 1;
	if (flag)
		board = new int*[ARRSIZE];
	for (int i = 0; i < ARRSIZE; i++)
	{
		if (flag)
			board[i] = new int[ARRSIZE];
		for (int j = 0; j < ARRSIZE; j++)
			board[i][j] = -1;
	}
}


void getPoint(const int chess_M, const int chess_N, const int* chess_top, int** chess_board, 
					 const int chess_noX, const int chess_noY, int lastX, int lastY, int& chess_ans_x, int& chess_ans_y)
{
	int new_havechess = 0;
	for (int i = 0; i < chess_M; i++)
		for (int j = 0; j < chess_N; j++)
			if (chess_board[i][j] > 0)
				new_havechess++;
	//	������
	
	if (new_havechess < havechess)
	{
#ifdef PRINT
		kengdieba << "wocao" << " " << M << " " << N << " " 
			<< chess_M << " " << chess_N << " "
			<< havechess << " " << new_havechess << " " << board << endl;
		if (board != NULL)
		for (int l = 0; l < ARRSIZE; l++)
		{
			for (int r = 0; r < ARRSIZE; r++)
				kengdieba << board[l][r] << " ";
			kengdieba << endl;
		}
		kengdieba << "+++++++++++++++############" << endl;
		
		for (int l = 0; l < ARRSIZE; l++)
		{
			for (int r = 0; r < ARRSIZE; r++)
				kengdieba << Board[l][r] << " ";
			kengdieba << endl;
		}
		kengdieba << "+++++++++++++++++!!!!!!!!!!!" << endl;
		for (int l = 0; l < chess_M; l++)
		{
			for (int r = 0; r < chess_N; r++)
				kengdieba << chess_board[l][r] << " ";
			kengdieba << endl;
		}
		kengdieba << "+++++++++++++++++++@@@@@@@@@@" << endl;
#endif
		initial_keng();
		while (size_pos)
		{
			tree[size_pos].round_sum = tree[size_pos].round_win = 0;
			tree[size_pos--].son.clear(); 	
		}
		tree[1].player = max_No;
		next_expand[1] = 0;
		size_pos = 1;
		root = 1;
	}
	//	�¾ֳ�ʼ��

	//VAR_C = MAXVAR;
	havechess = new_havechess;
	last_x = lastX;
	last_y = lastY;
	M = chess_M;
	N = chess_N;	
	noX = chess_noX;
	noY = chess_noY;
	for (int i = 0; i < N; i++)
		Top[i] = chess_top[i];
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			Board[i][j] = chess_board[i][j];
	//	������Ϣ
	
	start_time = clock();
	ucth();
	chess_ans_x = ans_x;
	chess_ans_y = ans_y;

#ifdef PRINT
	if (ans_x != chess_top[ans_y] - 1)
	{
		for (int i = 0; i < N; i++)
		file << chess_top[i] << " ";
		file << "	opop	" << endl;
		file << ans_x << " " << ans_y << " " << chess_top[ans_y] << endl;
		file << M << " " << N << endl;
		for (int l = 0; l < M; l++)
		{
			for (int r = 0; r < N; r++)
				file << board[l][r] << " ";
			file << endl;
		}
		file << "*******************" << endl;
	/*	for (int l = 0; l < ARRSIZE; l++)
		{
			for (int r = 0; r < ARRSIZE; r++)
				file << derta_sub[l][r] << " ";
			file << endl;
		}
		file << "+++++++++++++++++++" << endl;*/
	}
		kengdieba << "my choice:	" << ans_x << " " << ans_y << endl;
#endif
}

/*
int imitate(int player)
{
#ifdef PRINT
	//kengdieba << "imitate begin" << clock() << endl;
#endif
	int last = 0;
	for (int i = 0; i < N; i++)
		//if (derta_sub[top[i]][i] > 0)
	{
		int derta = 1;
		if ((i == noY) && (top[i] - 1 == noX))
			derta = 2;
		if (top[i] - derta < 0)
			continue;
		rand_seq[last++] = i;
		derta_seq[i] = derta;
	}
	random_shuffle(rand_seq, rand_seq + last);
	for (int ii = 0; ii < last; ii++)
	{
		int i = rand_seq[ii];
		int derta = derta_seq[i];
		if (player == max_No)
		{
//			board[top[i] - derta_sub[top[i]][i]][i] = 2;
//			if (machineWin(top[i] - derta_sub[top[i]][i], i, M, N, board))

			board[top[i] - derta][i] = 2;
			if (machineWin(top[i] - derta, i, M, N, board))		
				return max_No;
		}
		else
		{
			board[top[i] - derta][i] = 1;
			if (userWin(top[i] - derta, i, M, N, board))
				return min_No;
		}
		board[top[i] - derta][i] = 0;
	}
	//	��ʤ��չ�������ʤ��

	
	if (last == 0)
	{
#ifdef PRINT
	kengdieba << "bukengxue" << endl;
#endif
		return min_No;
	}


	int rol = rand_seq[0];
	top[rol] -= derta_seq[rol];		
//	top[rol] -= derta_sub[top[rol]][rol];	
	board[top[rol]][rol] = 2 - player;
#ifdef PRINT
	//kengdieba << "imitate end" << clock() << endl;
#endif
	return imitate(player ^ 1);
}
//	���ģ��
*/

void expand(int node_pos)
{
	int col = next_expand[node_pos];
	if (col == N - 1)
		next_expand[node_pos] = N;
	else
		for (int& i = ++next_expand[node_pos]; i < N; i++)
		{
			int derta = 1;
			if ((i == noY) && (top[i] - 1 == noX))
				derta = 2;
			if (top[i] - derta > -1)
				break;
		}
	//	������һ����չλ��

	tree[node_pos].son.push_back(++size_pos);
	tree[size_pos].player = tree[node_pos].player ^ 1; 
	tree[size_pos].col = col;
	top[col]--;
	if ((col == noY) && (top[col] == noX))
		top[col]--;
	board[top[col]][col] = 2 - tree[node_pos].player;
	//	����
	//	player �������߶��������ڵ�

	tree[size_pos].round_sum = 1;
	if (tree[node_pos].player == max_No)
	{
		if (machineWin(top[col], col, M, N, board))
		{
			tree[size_pos].round_win = 1;
			next_expand[size_pos] = next_expand[node_pos] = N;
			//	����Ѷ���������Ӳ�����չ
		}
	}
	else if (userWin(top[col], col, M, N, board))
	{
		tree[size_pos].round_win = -1;
		next_expand[size_pos] = next_expand[node_pos] = N;
		//	����Ѷ���������Ӳ�����չ
	}
	if (! tree[size_pos].round_win)
	{
		next_expand[size_pos] = get_first_expand(tree[size_pos].player);
		//	��һ������չλ��
		initial();
		if (imitate(tree[size_pos].player) == max_No)
			tree[size_pos].round_win = 1;
		else
			tree[size_pos].round_win = -1;
		//	������
	}

#ifdef PRINT
	/*kengdieba << "~~~~~~~~~~~~~~~" << endl;
	for (int i = 0; i < M; i++){
		for (int j = 0; j < N; j++)
			kengdieba << board[i][j] << " ";
		kengdieba << endl;}
	kengdieba << size_pos << " " << tree[node_pos].tcb << " " << tree[size_pos].col << " "
		<< tree[size_pos].round_win << " " << tree[size_pos].round_sum <<  " "
		<< tree[root].round_win << " " << tree[root].round_sum << endl;
	kengdieba << "~~~~~~~~~~~~~~~" << endl;*/
#endif
	rwin = tree[size_pos].round_win;
}
//	��չҶ�ӽڵ�


bool cmp(const int& A, const int& B)
{
	return tree[A].tcb > tree[B].tcb;
}

#ifdef PRINT
int maxvhh = 0;
#endif
void update()
{
#ifdef PRINT
	if (stack_top > maxvhh)
	{
		maxvhh = stack_top;
		kengdieba << maxvhh << "	uuuooo	" << endl;
	}
#endif
	while (stack_top)
	{
		int node_pos = dad_stack[stack_top--];
		tree[node_pos].round_win += rwin;
		tree[node_pos].round_sum++;	
		int son_size = tree[node_pos].son.size();
		int modify_c;
		if (tree[node_pos].player == max_No)
			modify_c = 1;
		else
			modify_c = -1;
		for (int i = 0; i < son_size; i++)
		{
			int son_pos = tree[node_pos].son[i];
			tree[son_pos].tcb = (double)(tree[son_pos].round_win * modify_c) / (double)tree[son_pos].round_sum + 
			//max(CONST_C + VAR_C, 0.3) * sqrt(2 * log(tree[node_pos].round_sum) / (double)tree[son_pos].round_sum);
			CONST_C * sqrt(2 * log(tree[node_pos].round_sum) / (double)tree[son_pos].round_sum);
		}
		sort(tree[node_pos].son.begin(), tree[node_pos].son.end(), cmp);		
	}
	//VAR_C -= 0.00002;
}
//	���½ڵ���Ϣ


void copyarr()
{
	for (int i = 0; i < M; i++)
		for (int j = 0; j < N; j++)
			board[i][j] = Board[i][j];
	for (int i = 0; i < N; i++)
		top[i] = Top[i];
}
//	���Ƶ�ǰ���


void iterate_imitate()
{
	copyarr();
	int node_pos;
	node_pos = choose();
	if (node_pos != -1)
		expand(node_pos);
	else
	{
		if (tree[dad_stack[stack_top]].player == min_No)
			rwin = 1;
		//	��Ӯ��
		else
			rwin = -1;
	}
	update();
}
//	��չģ��


int get_first_expand(int player)
{
	int re = N;
	for (int i = 0; i < N; i++)
		{
			int derta = 1;
			if ((i == noY) && (top[i] - 1 == noX))
				derta = 2;
			if (top[i] - derta > -1)
			{
				re = i;
				break;
			}
		}
	for (int i = 0; i < N; i++)
	{
		int derta = 1;
		if ((i == noY) && (top[i] - 1 == noX))
			derta = 2;
		if (top[i] - derta < 0)	
			continue;
		int must = -1;
		if (player == max_No)
		{
			board[top[i] - derta][i] = 2;
			if (machineWin(top[i] - derta, i, M, N, board))
				must = i;
		}
		else
		{
			board[top[i] - derta][i] = 1;
			if (userWin(top[i] - derta, i, M, N, board))
				must = i;
		}
		board[top[i] - derta][i] = 0;
		if (must != -1)
		{
			re = must;
			break;
		}
	}
	return re;
}
//	���ص�һ������չ��λ��


int choose()
{
	int node_pos = root, re = -1;
	stack_top = 0;
	//	��ʼ��
	while (re == -1)
	{
		dad_stack[++stack_top] = node_pos;
		if (next_expand[node_pos] < N)
			re = node_pos;
		else
		{
			int son_size = tree[node_pos].son.size();
			if (! son_size)
				return -1;
			for (int i = 1; i < son_size; i++)
				if (tree[tree[node_pos].son[i]].tcb < tree[tree[node_pos].son[0]].tcb)
				{
					son_size = i;
					break;
				}
			if (rand() < 30)
				son_size = tree[node_pos].son.size();
			node_pos = tree[node_pos].son[rand() % son_size];
			top[tree[node_pos].col]--;	
			if ((tree[node_pos].col == noY) && (top[tree[node_pos].col] == noX))
				top[tree[node_pos].col]--;	
			board[top[tree[node_pos].col]][tree[node_pos].col] = 1 + tree[node_pos].player;
			//	����Ӧ��������������
		}
	}
	return re;
}
//	����һ������չ�ڵ�


double ucth()
{
#ifdef PRINT
	//kengdieba << MAXVAR << " zhenima " << VAR_C << endl;
#endif
	copyarr();
	int have_next = 0;
	for (int i = 0; i < tree[root].son.size(); i++)
	{
		int son_pos = tree[root].son[i];
		if (tree[son_pos].col == last_y)
		{
			root = son_pos;
			have_next = 1;
			break;
		}
	}
	if (! have_next)
	{
		while (size_pos)
		{
			tree[size_pos].round_sum = tree[size_pos].round_win = 0;
			tree[size_pos--].son.clear(); 	
		}
		tree[1].player = max_No;
		next_expand[1] = get_first_expand(max_No);
		size_pos = 1;
		root = 1;
	}
#ifdef PRINT
	if (! have_next)
		kengdieba << "	nimabianimabia	" << endl;
#endif
	//	���¸�

	int iterate_last = 1000000000;	
	for (int iterate_time = 0; iterate_time < iterate_last; iterate_time++)
	{
		iterate_imitate();
		if (((double)(clock() - start_time) / 1000.0 > 4.5) || (size_pos > 8000000))
		{
			/*if (! keng_flag)
			{
				MAXVAR = iterate_time * 0.00002 - 2;
				keng_flag = 1;
			}*/
			break;
		}
#ifdef PRINT
		if (iterate_time % 10000 == 0)
		{
			kengdieba << "time past::::::::	" << iterate_time << " " << endl;
			kengdieba << (double)(clock() - start_time) / 1000.0 << endl;
			kengdieba << "++++++++++++++++++++" << endl;
		}
#endif
	}
	//	����

	int son_size = tree[root].son.size(), means = 0;
	for (int i = 0; i < son_size; i++)
	{
		int son_pos = tree[root].son[i];
		means += tree[son_pos].round_sum;
	}
	means /= son_size;
	int heheans[ARRSIZE], hehea = 1;
	for (int i = 0; i < son_size; i++)
	{
		int son_pos = tree[root].son[i];
		if (tree[son_pos].round_sum >= means)
		{
			heheans[0] = i;
			break;
		}
	}	
	for (int i = heheans[0] + 1; i < son_size; i++)
	{
		int son_pos = tree[root].son[i];
		if (tree[son_pos].tcb < tree[tree[root].son[heheans[0]]].tcb)
			break;
		else
		{
			if ((tree[son_pos].tcb != -1) && (tree[son_pos].round_sum >= means))
				heheans[hehea++] = i;
		}
	}
	int new_root = tree[root].son[heheans[rand() % hehea]];
	ans_y = tree[new_root].col;
	int derta = 1;
	if ((ans_y == noY) && (Top[ans_y] - 1 == noX))
		derta = 2;	
	ans_x = Top[ans_y] - derta;
	//	�ҳ�ʤ������
	

#ifdef PRINT
	kengdieba << "&&&&&&&&&&&&&&&&&&" << endl;
	kengdieba << tree[root].round_sum << " " << tree[root].round_win << endl;
	for (int i = 0; i < tree[root].son.size(); i++)
	{
		int son_pos = tree[root].son[i];
		kengdieba << i << "	" << tree[son_pos].col << "	" << son_pos << " " << tree[son_pos].tcb << "	" 
			<< tree[son_pos].round_sum << "	" << tree[son_pos].round_win 
			<< " " << (double)tree[son_pos].round_win / (double)tree[son_pos].round_sum << endl;
		for (int j = 0;j < tree[son_pos].son.size(); j++)
			kengdieba << tree[tree[son_pos].son[j]].col << " ";
		kengdieba << endl;
	}
	kengdieba << "&&&&&&&&&&&&&&&&&&" << endl;
	kengdieba << "uct	effort	" << endl;
	kengdieba << size_pos << "	" << tree[root].round_sum << "	" << tree[root].round_win << endl;
	kengdieba << ":::::::::::::::::::::" << endl;
#endif
	root = new_root;
	return tree[root].tcb;
}


int imitate(int player)
{
	int re = min_No;
	while (seq_pos < seq_point)
	{
		int rol = random_sequence[seq_pos++];
		int derta = 1;
		if ((rol == noY) && (top[rol] - 1 == noX))
			derta = 2;
		if (top[rol] - derta < 0)
			continue;
		top[rol] -= derta;
		if (player == max_No)
		{
			board[top[rol]][rol] = 2;
			if (machineWin(top[rol], rol, M, N, board))
				re = max_No;
			else
				re = imitate(player ^ 1);
		}
		else
		{
			board[top[rol]][rol] = 1;
			if (! userWin(top[rol], rol, M, N, board))
				re = imitate(player ^ 1);
		}
		break;
	}
	return re;
}
//	���ģ��


void initial()
{
	seq_point = 1;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < M; j++)
			random_sequence[seq_point++] = i;
	random_shuffle(random_sequence, random_sequence + seq_point);
	seq_pos = 1;
}
//	���ģ���������λ��
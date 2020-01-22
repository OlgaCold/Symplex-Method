#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M 20
#define N 20
#define IND 70

typedef struct{
    int m, n;   //number of column = n, rows = m 
    double matrix[M][N];    // matrix(MxN)

}Table;

bool equal (double a, double b);
void indent(int p_number);
bool read_table(Table *table, const char* p_filename);
bool print_table(Table *table, const char* type, int *basis_col);
void add_slack_variables(Table *table);
int find_pivot_column(Table *table, bool mode);
int find_pivot_row(Table *table, int pivot_col);
bool table_recalc(Table *table, int row, int col);
int find_basis_variable(Table *table, int column);
void print_solution(Table *table, const char* message);
bool symplex(Table *table, bool mode, int *basis_col, double *sol_vector);
int check_b_positive(Table *table);

const double epsilon = 1.0e-8;
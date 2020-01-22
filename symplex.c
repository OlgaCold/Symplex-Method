#include "symplex.h"

bool equal (double a, double b){
    return abs(a-b) < epsilon;
}

void indent(int p_number){    // отступ между таблицами
    for(int i = 0; i < p_number; i++){
        printf("-");
    }
    printf("\n");
}

bool read_table(Table *table, const char* p_filename){
    int input;
    FILE *fp;
 
    fp = fopen(p_filename, "r");
    if (fp == NULL){
        printf("Could not open file %s\n",p_filename);
        system("pause");
        return false;
    }
    input = fscanf(fp, "%d %d", &table->m, &table->n);
    if(table->m == 0 ||table->n == 0 || input == EOF || input == 0){ 
        printf("Cannot read dimensions\n"); 
        return false;
    }
    for(int i = 0; i < table->m; i++){
        for(int j = 0; j < table->n; j++){
            input = fscanf(fp, "%lf", &table->matrix[i][j]);
            if (input == EOF) {
                printf("Cannot read matrix\n"); 
                system("pause");
                return false;
            }
        }
    }
    fclose(fp);
    return 1;
}

bool print_table(Table *table, const char* type, int *basis_col){
    static int counter = 1; 
    indent(IND);
    printf("%d.  %s table: \n", counter++, type);
    indent(IND);
    printf("col:\tB[i]\t");
    for(int j = 1; j < table->n; j++){ printf("X%d\t", j); } 
    printf("\n");

    for(int j = 1; j < table->n; j++){
        int basis = find_basis_variable(table, j);
        if(basis != -1){
            basis_col[basis - 1] = j;
        }
    }

    for(int i = 0; i < table->m; i++){
        if(i == 0) printf("F:\t");
        //else printf("b%d:\t", i);
        else if(basis_col[i-1]){ printf("X%d:\t", basis_col[i-1]); }
        else { printf("\t"); }
        for(int j = 0; j < table->n; j++){
            if (equal((int)table->matrix[i][j], table->matrix[i][j]))
                printf("%d\t", (int)table->matrix[i][j]);
            else
                printf("%.2lf\t", table->matrix[i][j]);         
        }
        printf("\n");  
    }
    indent(IND);
}

void add_slack_variables(Table *table){
    for(int i = 1; i < table->m; i++){
        for(int j = 1; j < table->m; j++){
            table->matrix[i][j + table->n -1] = (i == j) ? 1: 0;
        }
    }
    table->n += table->m -1;
}

int find_pivot_column(Table *table, bool mode){
    if(!mode){//minimize
        int pivot_col = 1;
        double highest = table->matrix[0][pivot_col];
        for(int j = 1; j < table->n; j++){
            if(table->matrix[0][j] > highest){
                highest = table->matrix[0][j];
                pivot_col = j;
            }
        }
        
        if(highest < 0 || equal(highest,0)){
            printf("All coefficients in function are negative. It is optimal.\n");
            return -1;
        } else { 
            printf("Most positive column is X%d = %.2lf. It is pivot column - X%d.\n", pivot_col, highest, pivot_col);
            return pivot_col;

        }
    }else{//maximize
        int pivot_col = 1;
        double lowest = table->matrix[0][pivot_col];
        for(int j = 1; j < table->n; j++){
            if(table->matrix[0][j] < lowest){
                lowest = table->matrix[0][j];
                pivot_col = j;
            }
        }
        if(lowest > 0 ||  equal(lowest,0)){
            printf("All coefficients in function are positive. It is optimal.\n");
            return-1;
        } else {
            printf("Most negative column is X%d = %.2lf. It is pivot column - X%d.\n", pivot_col, lowest, pivot_col);
            return pivot_col;
        }
    }
}

int find_pivot_row(Table *table, int pivot_col){
    printf("To find pivot row to take out of basis,\nlets find ratios of free member of the row to an element of the pivot column:\n");
    int pivot_row = 0;//
    double min_ratio = -1;
    for(int i = 1; i < table->m; i++){
        double ratio = table->matrix[i][0]/table->matrix[i][pivot_col];
        printf("Matrix[%d][0]/matrix[%d][%d] = %.2lf,\n", i, i, pivot_col);
        if((ratio > 0  && ratio < min_ratio) || min_ratio < 0){
            min_ratio = ratio;
            pivot_row = i;
        } 
    }
    if (min_ratio < 0) return -1;
    printf("Min positive ratio = %.3lf in row %d. Found pivot row = %d.\n", min_ratio, pivot_row, pivot_row);
    return pivot_row;
}

int find_row_(Table *table, int pivot_col){
    printf("To find pivot row to take out of basis,\nlets find ratios of free member of the row to an element of the pivot column:\n");
    int pivot_row = 0;//
    double min_ratio = -1;
    for(int i = 1; i < table->m; i++){
        double ratio = table->matrix[i][0]/table->matrix[i][pivot_col];
        printf("Matrix[%d][0]/matrix[%d][%d] = %.2lf,\n", i, i, pivot_col);
        if((ratio > 0  && ratio < min_ratio) || min_ratio < 0){
            min_ratio = ratio;
            pivot_row = i;
        } 
    }
    if (min_ratio < 0) return -1;
    printf("Min positive ratio = %.3lf in row %d. Found pivot row = %d.\n", min_ratio, pivot_row, pivot_row);
    return pivot_row;
}

bool table_recalc(Table *table, int row, int col){
    double pivot_element = table->matrix[row][col];
    if(equal(pivot_element,0)){
        return false;
    }
    for(int j = 0; j < table->n; j++)
        table->matrix[row][j] /= pivot_element;
  
    for(int i = 0; i < table->m; i++) { // foreach remaining row i do
        double multiplier = table->matrix[i][col];
        if(i == row) continue;
        for(int j = 0; j < table->n; j++) {
            table->matrix[i][j] -= multiplier * table->matrix[row][j];
        }
    }
    return true;

}

int find_basis_variable(Table *table, int column){
    int basis_col = -1;
    for(int i = 1; i < table->m; i++){
        if(equal(table->matrix[i][column], 1)){
            if(basis_col == -1){
                basis_col = i; // found first '1', save this row number.
            } else {
                return -1; //found second '1', not an identity matrix.
            }
        } else if(!equal(table->matrix[i][column], 0)){
            return -1; // not an identity matrix column.
        }
    }
    return basis_col;
}

void print_solution(Table *table, const char* message, double *sol_vector){
    printf("%s is: F = %3.2lf, \n", message, table->matrix[0][0]);
    sol_vector[0] =  table->matrix[0][0];
    for(int j = 1; j < table->n; j++){
        int coef = find_basis_variable(table, j);
        if(coef != -1){
            printf("X%d = %.2lf,\t", j, table->matrix[coef][0]); 
            sol_vector[j] =  table->matrix[coef][0];
        }
        else{
            printf("X%d = 0,\t", j);
            sol_vector[j] = 0;
        }
    }
    printf("\n");
}

int check_b_positive(Table *table) {
    for(int i = 1; i < table->m; i++){
        if(table->matrix[i][0] < 0) return i;
    }
    return 0;
}

bool symplex(Table *table, bool mode, int *basis_col, double *sol_vector){
    int loop = 0;

    add_slack_variables(table);
    int n;
    while(!(n = check_b_positive(table))){
        
    }
    /*if(!check_b_positive(table)) {
        printf("There are zeros in free terms of equations column!!!");
        return false;
    }*/
    print_table(table, "Rewrited (with basis variables)", basis_col);
    
    while(++loop){

        int pivot_col, pivot_row;
        pivot_col = find_pivot_column(table, mode);
        if(pivot_col < 0){
            if(!check_b_positive(table)) {
                printf("There are negatives in free terms of equations column!!!");
                pivot_row = check_b_positive(table);
                pivot_row = 
                //return false;
            }
            print_solution(table, "The optimal solution", sol_vector);
            return true;
        }else{
            printf("The variable X%d will be entered in the basis.\n", pivot_col);
        }
        pivot_row = find_pivot_row(table, pivot_col);
        if(pivot_row < 0){
            printf("There aren't positive ratio -> function unbounded.\n");
            return false;
        }
        if(!table_recalc(table, pivot_row, pivot_col)){
            printf("Pivot element = 0!!!\n");
            return false;
        }
        print_table(table, "Recalculated (Gauss)", basis_col);

        if(loop > 20) {
            printf("Too many iterations > %d.\n", loop);
            return false;
            //break;
        }
    }
}

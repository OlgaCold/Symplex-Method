#include "symplex.h"
#include "branch.h"

int main(){

    bool mode;
    int temp;
    const char* filename = "ex.txt";  
    Table table;

    if(read_table(&table, filename) == false){
        return 0;
    }

    double sol_vect[table.n + table.m - 1];

    int basis_columns[table.m - 1];
    for(int i = 0; i < table.m - 1; i++){
        basis_columns[i] = 0;
    }

    printf("Please, choose mode: maximize(1) or minimize(0) function: \n");
    scanf("%d", &temp);
    mode = temp;
    printf("You choose %simization.\n", (mode == 0) ? "min":"max");

    print_table(&table, "Initial", basis_columns);
    
    symplex(&table, mode, basis_columns, sol_vect);
    system("pause");
    return 0;
} 
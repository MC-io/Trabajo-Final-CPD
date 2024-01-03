#include <stdio.h>
#include <stdlib.h>
#include <time.h>

FILE *fin;       
int nodes_count; 
int **matrix;  


// Algoritmo serial
void calculate_triangles()
{
    int triangles_count = 0;
    int **mult = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * nodes_count * nodes_count);
    if (mult == NULL) {
      printf("Error: malloc for mult failed.\n");
      exit(1);
    }
    int *ptr = (int*)(mult + nodes_count);
    for (int i = 0; i < nodes_count; i++) {
        mult[i] = (ptr + nodes_count * i);
        for (int j = 0; j < i; j++) {
            mult[i][j] = 0;    
            for (int k = 0; k < i; k++) {
                mult[i][j] += matrix[i][k] * matrix[k][j];
            }
            triangles_count += mult[i][j] * matrix[i][j];
        }
    }
    
    printf("El grafo contiene: %d triangulos\n", triangles_count);
    free(mult);    
}

void initialize_matrix()
{
    int i,j,fscanf_result;
    double w;
    matrix = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * nodes_count * nodes_count);
    if (matrix == NULL) {
      printf("Error: malloc for matrix failed.\n");
      exit(1);
    }
    int *ptr = (int*)(matrix + nodes_count);
    for(i = 0; i < nodes_count; i++) {
        matrix[i] = (ptr + nodes_count * i);        
        for (j = 0; j < nodes_count; j++) {
            matrix[i][j] = 0;            
        }
    }            
    fscanf_result = fscanf(fin, "%d", &i);
    while (fscanf_result != 1 || i != -1) {
        fscanf_result = fscanf(fin, "%d %lf \n", &j, &w);
        if (i != j && matrix[i][j] == 0) {
            if (i < j) {
                matrix[j][i] = 1;
            } else if (i > j) {
                matrix[i][j] = 1;
            }
        }                
        fscanf_result = fscanf(fin, "%d", &i);
    }
}


void syntax_message(char *compiled_name)
{
    printf("Correct syntax:\n");
    printf("%s <input-file> \n", compiled_name);
    printf("where: \n");
    printf("<input-file> is the file containing a generated graph by RandomGraph that the algorithm will use.\n");
}

int read_parameters(char **argv)
{
    char *input_filename = argv[1];
    if (input_filename == NULL) {
        printf("Input file parameter missing.\n");
        syntax_message(argv[0]);
        return 0;
    }
            
    fin = fopen(input_filename, "r");
    if (fin == NULL) {
        printf("Cannot open input file %s.\n", input_filename);
        return 0;        
    }

    printf("Contando los triangulos del grafo extraido de: %s\n", input_filename);
    return 1;
}

int main(int argc, char **argv)
{
    if (!read_parameters(argv)) {
        printf("Programa termina.\n");
        return -1;    
    }

    int fscanf_result = fscanf(fin, "%d \n", &nodes_count);
    if (fscanf_result != 1 || nodes_count > 0) {
        printf("Conteo de nodos: %d\n", nodes_count);
        printf("Algoritmo iniciado, por favor espere...\n");
        initialize_matrix();
        clock_t t1 = clock();
        calculate_triangles();            
        clock_t t2 = clock();
        printf("Algoritmo finalizado!\n");
        printf("Tiemoi usado: %f segs\n", ((float)t2 -t1) / CLOCKS_PER_SEC);    
    } else {
        printf("Archivo vacio.\n");        
    }
    fclose(fin);
    return 0;
}

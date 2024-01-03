
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

FILE *fin;                          // Archivo de entrada
int nodes_count;                    // Conteo de los nodos del grafo
int blocks_count;                   // Conteo de blqoues a procesar
int **matrix;                       // Matriz de adyacencia
int **block;                        // Bloque asignado
int **multiplication_blocks;        // Bloques recibidos para la multiplicacion de matrices
int **el_multiplication_blocks;     // Bloques recibidos para multiplicacion por elementos
int **qmap;                         // Mapeo de los bloques asignados a procesar
int q_i,q_j;                        // Posicion de procesos asignados
int rank;                           // Rank del proceso MPI
int size;                           // Tamano de proocesos MPI
int q;                              // Tamano de la tabla de bloques
int n;                              // Tamano del bloque
int mult_blocks_receive_counter;    // Conteo de bloques recibidos para multiplicacin de matrices.
int el_mult_blocks_receive_counter; // Conteo de bloques recibidos para multiplicacin de matrices.
int triangles_count;                // Triangulos encontrados por cada proceso
int nw;
// Algoritmo serial de conteo de triangulos
void calculate_triangles()
{
    int triangles_count = 0;
    int ** mult = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * nodes_count * nodes_count);
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

// Inicializar matriz de adyacencia leyendo el archivo
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

    for (int i = nw; i < nodes_count; i++)
    {
        matrix[i][0] = 1;
        matrix[0][i] = 1;
    }
}

// Mensaje de error por parametros
void syntax_message(char *compiled_name)
{
    printf("Correct syntax:\n");
    printf("%s <input-file> \n", compiled_name);
    printf("where: \n");
    printf("<input-file> is the file containing a generated graph by RandomGraph that the algorithm will use.\n");
}

// Lee el archivo de matriz de adyacencia y lo valida
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

    printf("Counting Triangles of Graph retrieved from input file: %s\n", input_filename);
    return 1;
}

// Reduce el conteo de traingulos locales al proceso master 
// Si todos los procesos fueron asignados con bloques se utiliza mpi reduce
// sino cada proceso manda su resultado local al master
void reduce_triangles_counts()
{
    int total_triangles=0;
    if (blocks_count==size) {
        MPI_Reduce(&triangles_count, &total_triangles, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);        
    } else {
        if (rank==0) {
            MPI_Status status;
            int process_triangles,i;
            for (i = 1; i < blocks_count; i++) {
                MPI_Recv(&process_triangles, 1, MPI_INT, i, 0, MPI_COMM_WORLD, &status);
                total_triangles += process_triangles;
            }
            total_triangles += triangles_count;
        } else {
            MPI_Send(&triangles_count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);        
        }        
    }
    if (rank == 0) {
        printf("Graph contains: %d triangles\n", total_triangles);    
    }    
}

void calculate_local_triangles_master()
{
    int i,j,k;
    for (i = 1; i < q; i++) {
        MPI_Send(&block[0][0], n*n, MPI_INT, qmap[i][0], 0, MPI_COMM_WORLD);
    }

    triangles_count = 0;
    int** mult = (int**)malloc(sizeof(int*) * nodes_count + sizeof(int) * n * n);
    if (mult == NULL) {
      printf("Error: malloc for mult failed.\n");
      exit(1);
    }
    int* ptrm = (int*)(mult + nodes_count);
    for (i=0; i<n; i++) {
        mult[i] = (ptrm + n * i);
        for (j = 0; j < n; j++) {
            mult[i][j] = 0;    
            for (k = 0; k < n; k++) {
                mult[i][j] += block[i][k] * block[k][j];    
            }
            triangles_count += mult[i][j] * block[i][j];
        }
    }
    
    free(mult);    
    free(matrix);
    free(block);    
    free(qmap);
}

void calculate_local_triangles_slave()
{
    int i,j,k,b;
    triangles_count = 0;
    int **mult = (int**)malloc(sizeof(int *) * n + sizeof(int) * n * n);
    if (mult == NULL) {
      printf("Error: malloc for mult failed.\n");
      exit(1);
    }
    int *mult_ptr = (int*)(mult + n);

    for (i = 0; i < n; i++) {
        mult[i] = (mult_ptr + n * i);
    }

    for (b = 0; b < mult_blocks_receive_counter; b++) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mult[i][j]=0;
                for (k = 0; k < n; k++) {
                    mult[i][j] += block[i][k] * multiplication_blocks[k][(b*n)+j];        
                }
                triangles_count += mult[i][j] * el_multiplication_blocks[i][(b*n)+j];
            }
        }
    }
    
    // Computando los bloques asignados faltantes para los processos en el Q map diagonal
    // ya que las matrices de los bloques de mult de matrices y elementos no son del mismo
    //tamano
    if (q_i == q_j) {
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                mult[i][j] = 0;    
                for (k = 0; k < n; k++) {
                    mult[i][j] += block[i][k] * block[k][j];    
                }
                triangles_count += mult[i][j] * block[i][j];
            }
        }    
    }
    
    free(mult);
    free(block);
    free(multiplication_blocks);
    free(el_multiplication_blocks);
    free(qmap);    
}

// En esta funcion, cada proceso esclavo recibe los bloques requeridos de otros procesos
// para realizar este calculo. Cada proceso debe encontrar que proceso tiene su bloque requerido
// basados en el Q Map. Cada proceso requiere bloques de la fila de su columna en el Q map
// para la multiplicacion de matrices, y los blqoues previos de su fila en el Q Map,
// para la multiplicacion por elementos
void receive_blocks_from_other_processes()
{
    int i,j,k;
    MPI_Status status;

    // Guardar memoria para recibir bloques requeridos
    multiplication_blocks = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * ((q_j+1)*n));
    if (multiplication_blocks == NULL) {
      printf("Error: malloc for multiplication_blocks failed.\n");
      exit(1);
    }
    int *multiplication_blocks_ptr = (int*)(multiplication_blocks + n);
    el_multiplication_blocks = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * ((q_j+1)*n));
    if (el_multiplication_blocks == NULL) {
      printf("Error: malloc for el_multiplication_blocks failed.\n");
      exit(1);
    }
    int *el_multiplication_blocks_ptr = (int*)(el_multiplication_blocks + n);
    int **buffer = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * n);
    if (buffer == NULL) {
      printf("Error: malloc for buffer failed.\n");
      exit(1);
    }
    int *buffer_ptr = (int*)(buffer + n);
    for (i = 0; i < n; i++) {
        multiplication_blocks[i] = (multiplication_blocks_ptr + ((q_j+1)*n) * i);
        el_multiplication_blocks[i] = (el_multiplication_blocks_ptr + ((q_j+1)*n) * i);
        buffer[i] = (buffer_ptr + n * i);    
    }

    // Recibe blqoues de la fila de proceso columna en el Q Map.
    // Los procesos de la columan 0 solo requieren los bloques del proceso master.
    mult_blocks_receive_counter=0;
    if (q_j == 0) {
        mult_blocks_receive_counter++;        
        MPI_Recv(*buffer, n*n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                multiplication_blocks[i][j] = buffer[i][j];
            }
        }
    } else {
        for (j = 0; j < q_j + 1; j++) {
            if (qmap[j][q_i] != rank) {
                mult_blocks_receive_counter++;    
                MPI_Recv(*buffer, n*n, MPI_INT, qmap[q_j][j], 0, MPI_COMM_WORLD, &status);
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        multiplication_blocks[i][(j*n)+k] = buffer[i][k];
                    }
                }
            }        
            
        }
    }
    
    //Enviar bloques asignados en la columna de proceso fila en el Q Map. 
    for (i = q_i; i < q; i++) {
        if (qmap[i][q_i] != rank) {
            MPI_Send(&block[0][0], n*n, MPI_INT, qmap[i][q_i], 0, MPI_COMM_WORLD);        
        }            
    }

    // Recibe bloques de los procesos previos en la fila de procesos en el Q Map.
    el_mult_blocks_receive_counter=0;
    for (j = 0; j < q_j; j++) {
        el_mult_blocks_receive_counter++;
        MPI_Recv(*buffer, n*n, MPI_INT, qmap[q_i][j], 1, MPI_COMM_WORLD, &status);
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                el_multiplication_blocks[i][(j*n)+k] = buffer[i][k];
            }
        }
    }

    // Copiar bloque asignado en la ultima posicion de el_multiplication_blocks.
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            el_multiplication_blocks[i][(el_mult_blocks_receive_counter*n)+j] = block[i][j];
        }
    }

    // Enviar bloque asignado a los siguientes procesos en la fila del proceso en el Q Map.
    for (j = q_j + 1; j < q_i + 1; j++) {
        MPI_Send(&block[0][0], n*n, MPI_INT, qmap[q_i][j], 1, MPI_COMM_WORLD);
    }
    free(buffer);
}

// Recibe el mapeo de Q desde el proceso maestro
void retrieve_q_mapping()
{
    int i,j;
    qmap = (int**)malloc(sizeof(int*) * q + sizeof(int) * q * q);
    if (qmap == NULL) {
      printf("Error: malloc for qmap failed.\n");
      exit(1);
    }
    int *qptr = (int*)(qmap + q);
    for (i = 0; i < q; i++) {
        qmap[i] = (qptr + q * i);
    }
    // Recibe q mapping del maestro
    MPI_Bcast(*qmap, q*q, MPI_INT, 0, MPI_COMM_WORLD);
    for (i = 0; i < q; i++) {
        for (j = 0; j < i+1; j++) {
            if (qmap[i][j] == rank) {
                q_i=i;
                q_j=j;
                i=q;
                j=q+1;
            }            
        }            
    }
}

// Recibe el bloque asignado del proceso
void receive_block()
{
    int i,j,k;
    MPI_Status status;
    // Recibe tamano del bloque para guardar memoria para recibir el bloque.
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    block = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * n);
    int* ptr = (int*)(block + n);
    for (i = 0; i < n; i++) {
        block[i] = (ptr + n * i);            
    }
    // Recibe bloque asignado por el maestro
    MPI_Recv(*block, n*n, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);    
}


// Usada por el proceso maestro, calcula el Q mapping y envia a cada proceso su bloque asignado
void scatter_blocks()
{
    int i,j,p,q_row,q_column,p_row,p_column,counter;
    
    qmap = (int**)malloc(sizeof(int*) * q + sizeof(int) * q * q);
    if (qmap == NULL) {
      printf("Error: malloc for qmap failed.\n");
      exit(1);
    }
    int *qptr = (int*)(qmap + q);
    for (i = 0; i < q; i++) {
        qmap[i] = (qptr + q * i);
    }

    // Hace broadcaste del tamano del bloque a los demas procesos
    n = nodes_count/q;
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);    
    block = (int**)malloc(sizeof(int*) * n + sizeof(int) * n * n);
    if (block == NULL) {
      printf("Error: malloc for block failed.\n");
      exit(1);
    }
    int *ptr = (int*)(block + n);
    for (i=0; i<n; i++) {
        block[i] = (ptr + n * i);
    }

    // Calcula el Q mapping, y envia a cada proceso su bloque asignado.
    // El proceso maestro colocara el Id de cada proceso sequencialmente de arriba hacia abajo.
    // Por ejemplo a un Q Map de 4x4 se vera asi:
    //             [0, 0, 0, 0]
    //             [1, 4, 0, 0]
    //             [2, 5, 7, 0]
    //             [3, 6, 8, 9]
    q_row = 1;
    q_column = 0;
    counter = 0;
    for (p = 1; p < size; p++) {
        if (counter < blocks_count-1) {
            p_row=0;
            p_column=0;
            for (i = (q_row*n); i < (q_row*n) + n; i++) {
                for (j = (q_column*n); j < (q_column*n) + n; j++) {
                    block[p_row][p_column] = matrix[i][j];
                    p_column++;
                }
                p_column=0;
                p_row++;        
            }
            MPI_Send(&block[0][0], n*n, MPI_INT, p, 0, MPI_COMM_WORLD);
            qmap[q_row][q_column] = p;
            if (q_row<q-1) {                    
                q_row++;
            } else {
                q_column++;
                q_row = q_column;                    
            }
            counter++;
        }                    
    }

    // Mater tiene el primer bloque de la matriz.
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            block[i][j] = matrix[i][j];
        }        
    }

    // Hace Broadcast del Q Mapping al resto de procesos.
    MPI_Bcast(&qmap[0][0], q*q, MPI_INT, 0, MPI_COMM_WORLD);
}

int main(int argc, char **argv)
{
    // Inicializacion de MPI
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    
    // Cada proceso calcula su valor de Q y la cantidad de bloques
    // de esa forma sabra si un blqoue se le sera asignado   
    q = (-1+sqrt(1+8*size))/2; // Resolviendo P = Q(Q + 1)/2 para hallar Q
    blocks_count=((q*(q-1))/2)+q;
    if (rank == 0) {
        // Chequeo de parametros
        if (!read_parameters(argv)) {
            printf("Program terminates.\n");
            MPI_Abort(MPI_COMM_WORLD, -1);    
        }
        int fscanf_result = fscanf(fin, "%d \n", &nodes_count); 
        if (fscanf_result != 1 || nodes_count > 0) {
            printf("Conteo de nodos: %d\n", nodes_count);
            // La matriz requiere estar particionada conformemente.
            if (nodes_count % q != 0) 
            {
                nw = nodes_count;
                nodes_count += q - (nodes_count % q);
            }
            if (nodes_count % q ==0) {
                printf("Algoritmo iniciado, por favor espere...\n");
                initialize_matrix();
                clock_t t1 = clock();
                // Cuando q = 1realiza el algoritmo serial
                if (q > 1) {
                    scatter_blocks(); 
                    calculate_local_triangles_master(); 
                    reduce_triangles_counts(); 
                } else {
                    calculate_triangles();
                }                                                
                clock_t t2 = clock();
                printf("Algoritmo Finalizado!\n");
                printf("Tiempo usado: %f segs\n", ((float)t2 -t1) / CLOCKS_PER_SEC);
                fclose(fin);        
            } else {
                printf("El grafo no puede ser particionado\n");
                printf("Utiliza un numero de procesos diferente.\n");
                fclose(fin);
                MPI_Abort(MPI_COMM_WORLD, -1);
            }
        } else {
            printf("Archivo vacio.\n");
            fclose(fin);
            MPI_Abort(MPI_COMM_WORLD, -1);
        }    
    } else if (q > 1 && rank < blocks_count) {
        receive_block(); 
        retrieve_q_mapping();
        receive_blocks_from_other_processes(); 
        calculate_local_triangles_slave(mult_blocks_receive_counter); 
        reduce_triangles_counts();
    }
    MPI_Finalize();
}

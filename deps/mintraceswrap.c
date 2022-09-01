#include <traces.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

TracesOptions defaultoptions_traces() {
	DEFAULTOPTIONS_TRACES(options);
	return options;
}


static int * generator_array;
static size_t array_len;
static size_t used_len;
static bool * visited;
static size_t visited_len;

void recordautomproc(int count, int* perm, int n){
	for(int i = 0; i < n; i++){
		visited[i] = false;
	}
	for(int i = 0; i < n; i++){
		if (perm[i] == i || visited[i]){ // perm[i] == i, skip;  perm[i] < i, already recorded
			continue;
		}
		else
		{
			generator_array[used_len] = i;
			visited[i] = true;

			used_len++;
			if(used_len >= array_len){
				DYNREALLOC(int, generator_array, array_len, array_len + 1000, "recordautomproc: add generators");
			}

			int j = perm[i];
			while(j != i){
				generator_array[used_len] = j;
				visited[j] = true;
				
				
				used_len++;
				if(used_len >= array_len){
					DYNREALLOC(int, generator_array, array_len, array_len + 1000, "recordautomproc: add generators");
				}

				j = perm[j];
			}
			generator_array[used_len] =-1;

			used_len++;
			if(used_len >= array_len){
				DYNREALLOC(int, generator_array, array_len, array_len + 1000, "recordautomproc: add generators");
			}
		
		}
	}
	generator_array[used_len] = -2;

	used_len++;
	if(used_len >= array_len){
		DYNREALLOC(int, generator_array, array_len, array_len + 1000, "recordautomproc: add generators");
	}
}



int * Traces_With_Automs(sparsegraph* g, size_t * generator_len, int* labelling, int* partition, int* orbits, TracesOptions* options, TracesStats* stats, sparsegraph* outgraph){	
	array_len = 0;
	visited_len = 0;
	generator_array = NULL;
	visited = NULL;
	DYNALLOC1(bool, visited, visited_len, visited_len + g->nv, "Traces_With_Automs: allocate visited");
	DYNALLOC1(int, generator_array, array_len, array_len + 1000, "Traces_With_Automs: allocate generators");
	used_len = 0;
	options->userautomproc = recordautomproc;
	Traces(g, labelling, partition, orbits,  options, stats, outgraph);
	*generator_len = used_len;
	return generator_array;
}


void Traces_With_Automs_Free(){
	DYNFREE(generator_array, array_len);
	DYNFREE(visited, visited_len);
}


void Traces_With_Automs_DEBUG(sparsegraph* g, int* labelling, int* partition, int* orbits, TracesOptions* options, TracesStats* stats, sparsegraph* outgraph){
	options->writeautoms = TRUE;
	options->cartesian = FALSE;
	options->linelength = 0;
	FILE *stream;
	char *buf = NULL;
	size_t len = 0;
	stream = open_memstream(&buf, &len);
	if (stream == NULL)
		fprintf(ERRFILE, "cannot open stream for writing automorphisms\n");
	options->outfile = stream;
	Traces(g, labelling, partition, orbits, options, stats, outgraph);
	fflush(stream);
	fclose(stream);
	printf("%s %zu\n", buf, len);
	free(buf);	
}

/**
 * Copyright (C) 2012-2015 Universidad Simón Bolívar
 *
 * Copying: GNU GENERAL PUBLIC LICENSE Version 2
 * @author Guillermo Palma <gpalma@ldc.usb.ve>
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <time.h>

#include "types.h"
#include "memory.h"
#include "util.h"
#include "input.h"
#include "semEP.h"

#define MIN_ARG   7

struct global_args {
     char *matrix_left_filename;
     char *matrix_right_filename;
     char *left_filename;
     char *right_filename;
     char *graph_filename;
     double threshold_left;
     double threshold_right;
};

static struct global_args g_args;
static const char *optString = "l:r:";

/*********************************
 **  Parse Arguments
 *********************************/

static void display_usage(void)
{
     fatal("Incorrect arguments \n\tsemEP  [-l left threshold] [-r right threshold] <left matrix>  <left vertices> <right matrix> <right vertices> <bipartite graph>\n");
}

static void get_name(const char *filename, int n, char *instance)
{
     char buf[n];
     char *aux = NULL;

     strcpy(buf,filename);
     aux = strtok(buf,"/");
     do {
	  strcpy(instance, aux);
     } while((aux = strtok(NULL,"/")) != NULL);

     strcpy(buf,instance);
     aux = strtok(buf,".");
     strcpy(instance, aux);
}

static void initialize_arguments(void)
{
     g_args.matrix_left_filename = NULL;
     g_args.matrix_right_filename = NULL;
     g_args.left_filename = NULL;
     g_args.right_filename = NULL;
     g_args.graph_filename = NULL;
     g_args.threshold_left = 0.0;
     g_args.threshold_right = 0.0;
}

static void print_args(void)
{
     printf("\n**********************************************\n");
     printf("Parameters:\n");
     printf("Left matrix file name: %s\n", g_args.matrix_left_filename);
     printf("Left vertices file name: %s\n", g_args.left_filename);
     printf("Right matrix file name: %s\n", g_args.matrix_right_filename);
     printf("Right vertices file name: %s\n", g_args.right_filename);
     printf("Left threshold: %.3f\n", g_args.threshold_left);
     printf("Right threshold: %.3f\n", g_args.threshold_right);
     printf("Graph bipartite file name: %s\n", g_args.graph_filename);
     printf("************************************************\n");
}

static void parse_args(int argc, char **argv)
{
     int i, opt;

     if (argc != (MIN_ARG+3))
	  display_usage();
     
     initialize_arguments();
     opt = getopt(argc, argv, optString);
     while(opt != -1) {
	  switch(opt) {
	  case 'r':
	       g_args.threshold_right = strtod(optarg, (char **)NULL);
	       break;
	  case 'l':
	       g_args.threshold_left = strtod(optarg, (char **)NULL);
	       break;
	  case '?':
	       display_usage();
	       break;
	  default:
	       /* You won't actually get here. */
	       fatal("?? getopt returned character code 0%o ??\n", opt);
	  }
	  opt = getopt(argc, argv, optString);
     }
     i = optind;
     g_args.matrix_left_filename = argv[i++];
     g_args.left_filename = argv[i++];
     g_args.matrix_right_filename = argv[i++];
     g_args.right_filename = argv[i++];
     g_args.graph_filename = argv[i];
}

/*********************************
 *********************************
 **
 **       Main section
 **
 *********************************
 **********************************/

int main(int argc, char **argv)
{
     int len;
     clock_t ti, tf;
     static char *name;
     struct input_data in;
     double density;

     ti = clock();
     parse_args(argc, argv);
     print_args();
     len = strlen(g_args.graph_filename) + 1;
     name = xcalloc(len, 1);
     get_name(g_args.graph_filename, len, name);
     printf("\n**** Go semEP! **** \n");
     in = get_input_data(g_args.matrix_left_filename, g_args.matrix_right_filename,
			 g_args.left_filename, g_args.right_filename,
			 g_args.graph_filename);
     density = semEP_solver(&in.left_matrix, &in.right_matrix, &in.left_terms, &in.right_terms,
			    &in.td, &in.bpgraph, g_args.threshold_left, g_args.threshold_right, name);
     printf("Average density of the partitions: %.4f \n", density);
     printf("*** semEP finished ***\n");
     tf = clock();
     printf("\nTotal time %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
     free(name);
     free_input_data(&in);

     return 0;
}

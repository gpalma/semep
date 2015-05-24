/**
 * Copyright (C) 2012, 2013, 2014 Universidad Simón Bolívar
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
#include "graph.h"
#include "util.h"
#include "input.h"
#include "semEP.h"

#define MIN_ARG   4

struct global_args {
     char *graph_filename;
     char *desc_filename;
     char *annt1_filename;
     char *annt2_filename;
     double threshold1;
     double threshold2;
     double threshold3;
     bool prediction;
     bool matrix;
     bool description;
     enum measure d;
};

static struct global_args g_args;
static const char *optString = "dpmt:r:s:e:";

/*********************************
 **  Parse Arguments
 *********************************/

static void display_usage(void)
{
     fatal("Incorrect arguments \n\tsemEP [-e tax|str|ps] [-m] [-p] [-d] [-t threshold A1] [-r threshold A2] [-s threshold BT] <source sim> <descrp> <A1> <A2>\n");
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
     g_args.graph_filename = NULL;
     g_args.desc_filename = NULL;
     g_args.annt1_filename = NULL;
     g_args.annt2_filename = NULL;
     g_args.threshold1 = 0.25;
     g_args.threshold2 = 0.25;
     g_args.threshold3 = 0.25;
     g_args.prediction = false;
     g_args.matrix = false;
     g_args.description = false;
     g_args.d = DTAX;
}

static void print_args(void)
{
     printf("\n*************************************\n");
     printf("Parameters:\n");
     if (g_args.matrix)
	  printf("Matrix: %s\n", g_args.graph_filename);
     else
	  printf("Graph: %s\n", g_args.graph_filename);
     printf("Terms description: %s\n", g_args.desc_filename);
     printf("Annotations of Entity 1: %s\n", g_args.annt1_filename);
     printf("Annotations of Entity 2: %s\n", g_args.annt2_filename);
     printf("Threshold E1: %.3f\n", g_args.threshold1);
     printf("Threshold E2: %.3f\n", g_args.threshold2);
     printf("Threshold BT: %.3f\n", g_args.threshold3);
     printf("Get predicted links: %s\n", g_args.prediction ? "true" : "false");
     printf("Matrix input: %s\n", g_args.matrix ? "true" : "false");
     printf("Get the description of the annotations: %s\n", g_args.description ? "true" : "false");
     if (g_args.d == DTAX) {
	  printf("Measure: d_tax\n");
     } else if (g_args.d == DSTR) {
	  printf("Measure: d^str_tax\n");
     } else if (g_args.d == DPS) {
	  printf("Measure: d_ps\n");
     } else {
	  fatal("Unknown measure");
     }
     printf("*************************************\n");
}

static void parse_args(int argc, char **argv)
{
     int i, opt;

     initialize_arguments();
     opt = getopt(argc, argv, optString);
     while(opt != -1) {
	  switch(opt) {
	  case 'd':
	       g_args.description = true;
	       break;
	  case 'p':
	       g_args.prediction = true;
	       break;
	  case 'm':
	       g_args.matrix = true;
	       break;
	  case 'r':
	       g_args.threshold2 = strtod(optarg, (char **)NULL);
	       break;
	  case 's':
	       g_args.threshold3 = strtod(optarg, (char **)NULL);
	       break;
	  case 't':
	       g_args.threshold1 = strtod(optarg, (char **)NULL);
	       break;
	  case 'e':
	       if (strcmp(optarg, "tax") == 0) {
		    g_args.d = DTAX;
	       } else if (strcmp(optarg, "str") == 0) {
		    g_args.d = DSTR;
	       } else if (strcmp(optarg, "ps") == 0) {
		    g_args.d = DPS;
	       } else {
		    display_usage();
	       }
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
     if ((argc - optind) != MIN_ARG)
	  display_usage();
     i = optind;
     g_args.graph_filename = argv[i++];
     g_args.desc_filename = argv[i++];
     g_args.annt1_filename = argv[i++];
     g_args.annt2_filename = argv[i];
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
     static char *name1, *name2;
     struct input_data in;
     double sim;

     ti = clock();
     parse_args(argc, argv);
     print_args();

     /* get the names of the concepts */
     len = strlen(g_args.annt1_filename) + 1;
     name1 = xcalloc(len, 1);
     get_name(g_args.annt1_filename, len, name1);

     len = strlen(g_args.annt2_filename) + 1;
     name2 = xcalloc(len, 1);
     get_name(g_args.annt2_filename, len, name2);

     /* start solver */
     printf("\n**** semEP Begins ***\n");
     
     in = get_input_data(g_args.graph_filename,
			 g_args.desc_filename,
			 g_args.annt1_filename,
			 g_args.annt2_filename,
			 g_args.matrix, g_args.description);

     sim = annotation_partition(in.object, in.n, &in.anntt1, &in.anntt2,
				g_args.threshold1, g_args.threshold2, g_args.threshold3, 
				name1, name2, in.descriptions,
				g_args.prediction, g_args.matrix, g_args.d);
     printf("Average similarity of the partitions: %.4f \n", sim);
     printf("*** semEP Finished ***\n");
     tf = clock();
     printf("\nTotal Time %.3f secs\n", (double)(tf-ti)/CLOCKS_PER_SEC);
     free(name1);
     free(name2);
     free_input_data(&in, g_args.matrix);
}

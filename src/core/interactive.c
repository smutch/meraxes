#include "meraxes.h"

static int prompt_char(const char* message)
{
    int cont;

    if (run_globals.mpi_rank == 0) {
        char c;
        int tmp;

        printf("\n%s (y/n): ", message);
        fflush(stdout);

        c = (char)getchar();
        do {
            tmp = getchar();
        } while (tmp != '\n' && tmp != EOF);

        switch (c) {
        case 'y':
        case 'Y': {
            cont = 1;
            break;
        }

        case 'n':
        case 'N': {
            cont = 0;
            break;
        }

        default: {
            printf("\nI do not understand your input...\n");
            cont = -1;
            break;
        }
        }
    }

    MPI_Bcast(&cont, 1, MPI_INT, 0, run_globals.mpi_comm);
    return cont;
}

void continue_prompt(char* param_file)
{
    int rerun = -1;

    fflush(stderr);
    fflush(stdout);

    while (rerun < 0)
        rerun = prompt_char("Reread input file and rerun model?");

    MPI_Barrier(run_globals.mpi_comm);

    if (rerun) {
        if (run_globals.params.Flag_PatchyReion)
            init_reion_grids();

        read_parameter_file(param_file, 1);
    } else {
        run_globals.params.FlagInteractive = 0;
        MPI_Bcast(&(run_globals.params.FlagInteractive), 1, MPI_INT, 0, run_globals.mpi_comm);
    }

    return;
}

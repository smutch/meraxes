#include "parse_paramfile.h"
#include "core/misc_tools.h"
#include <assert.h>
#include <ctype.h>
#include <meraxes.h>
#include <regex.h>

static int compile_regex(regex_t* reg, const char* regex_text)
{
  int status = regcomp(reg, regex_text, REG_EXTENDED | REG_NEWLINE);

  if (status != 0) {
    char error_message[256];
    regerror(status, reg, error_message, 256);
    mlog_error("Regex error compiling '%s': %s", MLOG_MESG, regex_text, error_message);
    return 1;
  }
  return 0;
}

static int match_regex(regex_t* reg, const char* match_str, entry_t* entry)
{
  // "p" is a pointer into the string which points to the end of the previous match.
  const char* p = match_str;
  // "n_matches" is the maximum number of matches allowed.
  const int n_matches = 3;
  // "match" contains the matches found.
  regmatch_t match[n_matches];

  int start;
  int finish;

  while (1) {
    int nomatch = regexec(reg, p, n_matches, match, 0);
    if (nomatch)
      return nomatch;
    for (int ii = 0; ii < n_matches; ii++) {
      if (match[ii].rm_so == -1)
        break;
      start = (int)(match[ii].rm_so + (p - match_str));
      finish = (int)(match[ii].rm_eo + (p - match_str));
      if (ii == 1)
        sprintf(entry->key, "%.*s", (finish - start), match_str + start);
      if (ii == 2) {
        sprintf(entry->value, "%.*s", (finish - start), match_str + start);

        // Trim trailing whitespace
        for (int jj = finish - start - 1; jj >= 0; jj--) {
          if (isspace(entry->value[jj]))
            entry->value[jj] = '\0';
          else
            break;
        }
      }
    }
    p += match[0].rm_eo;
  }
}

int parse_paramfile(char* fname, entry_t entry[PARAM_MAX_ENTRIES])
{
  char buffer[PARAM_MAX_LINE_LEN];
  FILE* fin;
  int level_change;
  int last_level;
  int counter;
  regex_t reg;
  const char* regex_text = "[[:space:]]*([^[:space:]^#]+)[[:space:]]*:[[:space:]]*([^#]*)";

  if ((fin = fopen(fname, "r")) == NULL) {
    mlog_error("Failed to open parameter file %s\t...ABORTING", fname);
    ABORT(EXIT_FAILURE);
  }
  compile_regex(&reg, regex_text);

  for (int ii = 0; ii < PARAM_MAX_ENTRIES; ii++)
    entry[ii].level = 0;

  counter = 0;
  level_change = 0;
  last_level = 0;
  while (fgets(buffer, PARAM_MAX_LINE_LEN, fin) != NULL) {
    entry[counter].key[0] = '\0';
    entry[counter].value[0] = '\0';
    entry[counter].level = last_level;
    entry[counter].level += level_change;
    last_level = entry[counter].level;

    level_change = 0;
    for (int ii = 0; buffer[ii]; ii++) {
      level_change += (buffer[ii] == '{');
      level_change -= (buffer[ii] == '}');
    }
    buffer[strcspn(buffer, "{")] = '\0';
    buffer[strcspn(buffer, "}")] = '\0';
    buffer[strcspn(buffer, "#")] = '\0';
    buffer[strcspn(buffer, "\n")] = '\0';

    if (buffer[0] != '\0') {
      match_regex(&reg, buffer, &(entry[counter]));
      counter++;
    }
  }

  regfree(&reg);
  fclose(fin);

  return counter;
}

static inline int chars_to_snap(const char* string, int maxsnaps)
{
  int val = atoi(string);
  if (val < 0)
    val += maxsnaps;
  return val;
}

int parse_slices(const char* string, const int arrmax, int** indices)
{
  char sep[] = ", ";
  char* context_outer = NULL;
  char parse_string[STRLEN];

  sprintf(parse_string, "%s", string);
  char* p = strtok_r(parse_string, sep, &context_outer);

  // count the number of requested snapshots
  int count = 0;
  while (p != NULL) {
    // check to see if this token is a slice
    if (strchr(p, ':')) {
      char* context_inner = NULL;
      char* sp = NULL;
      int slice_count = 0;
      int slice[3] = { 0, arrmax, 1 };

      if (p[0] == ':')
        slice_count++;
      for (sp = strtok_r(p, ":", &context_inner); sp; sp = strtok_r(NULL, ":", &context_inner), slice_count++) {
        slice[slice_count] = chars_to_snap(sp, arrmax);
      }

      count += (slice[1] - slice[0]) / slice[2] - 1;
    }
    count++;
    p = strtok_r(NULL, sep, &context_outer);
  }

  // allocate the ListOutputSnaps array
  assert(count > 0);
  *indices = calloc(count, sizeof(int));

  // reset the token string and parse this time
  sprintf(parse_string, "%s", string);
  p = strtok_r(parse_string, sep, &context_outer);

  count = 0;
  while (p != NULL) {
    // check to see if this token is a slice
    if (strchr(p, ':')) {
      char* context_inner = NULL;
      char* sp = NULL;
      int slice_count = 0;
      int slice[3] = { 0, arrmax, 1 };

      if (p[0] == ':')
        slice_count++;
      for (sp = strtok_r(p, ":", &context_inner); sp; sp = strtok_r(NULL, ":", &context_inner), slice_count++)
        slice[slice_count] = chars_to_snap(sp, arrmax);

      for (int ii = slice[0]; ii < slice[1]; ii++, count++)
        (*indices)[count] = ii;

    } else {
      (*indices)[count++] = chars_to_snap(p, arrmax);
    }
    p = strtok_r(NULL, sep, &context_outer);
  }

  // sort the list and remove any duplicates
  qsort(*indices, count, sizeof(int), compare_ints);
  int jj = 0;
  for (int ii = 1; ii < count; ii++) {
    if ((*indices)[ii] != (*indices)[jj])
      (*indices)[++jj] = (*indices)[ii];
  }
  count = jj + 1;

  return count;
}

void parse_output_snaps(const char* string)
{
  if (run_globals.mpi_rank == 0) {
    int* snaplist;
    int count = parse_slices(string, run_globals.params.SnaplistLength, &snaplist);

    run_globals.NOutputSnaps = count;
    run_globals.LastOutputSnap = snaplist[count - 1];

#ifdef CALC_MAGS
    if (count != MAGS_N_SNAPS) {
      mlog_error("Number of entries in output snaplist does not match MAGS_N_SNAPS!");
      ABORT(EXIT_FAILURE);
    }
#endif

    run_globals.ListOutputSnaps = NULL;
    run_globals.ListOutputSnaps = malloc(sizeof(int) * count);
    assert(run_globals.ListOutputSnaps != NULL);

    memcpy(run_globals.ListOutputSnaps, snaplist, sizeof(int) * count);
    free(snaplist);

#ifdef DEBUG
    mlog("Parsed snaplist = [", MLOG_MESG);
    for (int ii = 0; ii < count; ++ii) {
      mlog(" %d", MLOG_CONT, run_globals.ListOutputSnaps[ii]);
    }
    mlog(" ]", MLOG_CONT);
#endif
  }

  // broadcast the data to all other ranks
  MPI_Bcast(&(run_globals.NOutputSnaps), 1, MPI_INT, 0, run_globals.mpi_comm);

  if (run_globals.mpi_rank > 0)
    run_globals.ListOutputSnaps = malloc(sizeof(int) * run_globals.NOutputSnaps);

  MPI_Bcast(run_globals.ListOutputSnaps, run_globals.NOutputSnaps, MPI_INT, 0, run_globals.mpi_comm);
  MPI_Bcast(&(run_globals.LastOutputSnap), 1, MPI_INT, 0, run_globals.mpi_comm);
}

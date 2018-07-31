#include <ctype.h>
#include <meraxes.h>
#include <parse_paramfile.h>
#include <regex.h>

static int compile_regex(regex_t* reg, const char* regex_text)
{
    int status = regcomp(reg, regex_text, REG_EXTENDED | REG_NEWLINE);

    if (status != 0) {
        char error_message[256];
        regerror(status, reg, error_message, 256);
        mlog_error("Regex error compiling '%s': %s",
            MLOG_MESG, regex_text, error_message);
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
                sprintf(entry->key, "%.*s", (finish - start),
                    match_str + start);
            if (ii == 2) {
                sprintf(entry->value, "%.*s", (finish - start),
                    match_str + start);

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

        // DEBUG
        // fprintf(stderr, "buffer = %s :: level_change = %d :: level = %d\n", buffer, level_change, entry[counter].level);

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

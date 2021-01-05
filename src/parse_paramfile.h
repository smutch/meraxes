#define PARAM_MAX_ENTRIES 200
#define PARAM_MAX_LINE_LEN 512
#define PARAM_TYPE_INT 801
#define PARAM_TYPE_FLOAT 802
#define PARAM_TYPE_DOUBLE 803
#define PARAM_TYPE_STRING 804
#define PARAM_TYPE_LONGLONG 805
#define PARAM_TYPE_UNUSED 809

typedef struct entry_t entry_t;
struct entry_t {
    int level;
    char key[64];
    char value[256];
};

int parse_paramfile(char* fname, entry_t entry[PARAM_MAX_ENTRIES]);

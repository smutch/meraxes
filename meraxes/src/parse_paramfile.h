#define PARAM_MAX_ENTRIES        100
#define PARAM_MAX_LINE_LEN       512
#define PARAM_TYPE_INT           801
#define PARAM_TYPE_FLOAT         802
#define PARAM_TYPE_DOUBLE        803
#define PARAM_TYPE_STRING        804

typedef struct entry_t entry_t;
struct entry_t
{
  int level;
  char key[64];
  char value[256];
};

int parse_paramfile(char *fname, entry_t entry[PARAM_MAX_ENTRIES]);

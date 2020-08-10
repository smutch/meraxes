#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void cn_quote()
{
  const size_t max_quote_size = 512;
  char quote_list[28][max_quote_size];
  int ind;

  srand((unsigned int)time(NULL));

  strncpy(quote_list[0],
          "CHUCK_NORRIS_EXCEPTION: When Chuck Norris throws exceptions, it’s across the room.",
          max_quote_size);
  strncpy(quote_list[1],
          "CHUCK_NORRIS_EXCEPTION: All arrays Chuck Norris declares are of infinite size, because Chuck Norris knows "
          "no bounds.",
          max_quote_size);
  strncpy(
    quote_list[2],
    "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t have disk latency because the hard drive knows to hurry the hell up.",
    max_quote_size);
  strncpy(quote_list[3], "CHUCK_NORRIS_EXCEPTION: Chuck Norris writes code that optimizes itself.", max_quote_size);
  strncpy(quote_list[4],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris can’t test for equality because he has no equal.",
          max_quote_size);
  strncpy(quote_list[5],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t need garbage collection because he doesn’t call .Dispose(), he "
          "calls .DropKick().",
          max_quote_size);
  strncpy(quote_list[6], "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s first program was kill -9.", max_quote_size);
  strncpy(quote_list[7], "CHUCK_NORRIS_EXCEPTION: Chuck Norris burst the dot com bubble.", max_quote_size);
  strncpy(quote_list[8],
          "CHUCK_NORRIS_EXCEPTION: All browsers support the hex definitions #chuck and #norris for the colors black "
          "and blue.",
          max_quote_size);
  strncpy(quote_list[9],
          "CHUCK_NORRIS_EXCEPTION: MySpace actually isn’t your space, it’s Chuck’s (he just lets you use it).",
          max_quote_size);
  strncpy(quote_list[10],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris can write infinite recursion functions... and have them return.",
          max_quote_size);
  strncpy(
    quote_list[11], "CHUCK_NORRIS_EXCEPTION: Chuck Norris can solve the Towers of Hanoi in one move.", max_quote_size);
  strncpy(quote_list[12], "CHUCK_NORRIS_EXCEPTION: The only pattern Chuck Norris knows is God Object.", max_quote_size);
  strncpy(quote_list[13], "CHUCK_NORRIS_EXCEPTION: Chuck Norris finished World of Warcraft.", max_quote_size);
  strncpy(quote_list[14],
          "CHUCK_NORRIS_EXCEPTION: Project managers never ask Chuck Norris for estimations…ever.",
          max_quote_size);
  strncpy(quote_list[15],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t use web standards as the web will conform to him.",
          max_quote_size);
  strncpy(quote_list[16],
          "CHUCK_NORRIS_EXCEPTION: “It works on my machine” always holds true for Chuck Norris.",
          max_quote_size);
  strncpy(quote_list[17],
          "CHUCK_NORRIS_EXCEPTION: Whiteboards are white because Chuck Norris scared them that way.",
          max_quote_size);
  strncpy(quote_list[18],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t do Burn Down charts, he does Smack Down charts.",
          max_quote_size);
  strncpy(quote_list[19], "CHUCK_NORRIS_EXCEPTION: Chuck Norris can delete the Recycling Bin.", max_quote_size);
  strncpy(quote_list[20], "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s beard can type 140 wpm.", max_quote_size);
  strncpy(quote_list[21],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris can unit test entire applications with a single assert.",
          max_quote_size);
  strncpy(quote_list[22],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t bug hunt as that signifies a probability of failure, he goes "
          "bug killing.",
          max_quote_size);
  strncpy(
    quote_list[23],
    "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s keyboard doesn’t have a Ctrl key because nothing controls Chuck Norris.",
    max_quote_size);
  strncpy(quote_list[24],
          "CHUCK_NORRIS_EXCEPTION: When Chuck Norris is web surfing websites get the message 'Warning: Internet "
          "Explorer has deemed this user to be malicious or dangerous. Proceed?''.",
          max_quote_size);
  strncpy(quote_list[25],
          "CHUCK_NORRIS_EXCEPTION: You cannot throw Chuck Norris... Chuck Norris, throws YOU!",
          max_quote_size);
  strncpy(quote_list[26],
          "CHUCK_NORRIS_EXCEPTION: Chuck Norris wrote Hello World once... it was called Unix.",
          max_quote_size);
  strncpy(quote_list[27], "CHUCK_NORRIS_EXCEPTION: Chuck Norris will asynchronously KICK YOUR ASS.", max_quote_size);

  ind = (int)((double)rand() * (28.0 + 1.0) / RAND_MAX);

  fprintf(stderr, "%s\n", quote_list[ind]);
  fflush(stderr);
}

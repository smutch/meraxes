#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void cn_quote()
{
  const size_t max_quote_size = 512;
  char quote_list[104][max_quote_size];
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
  strncpy(quote_list[28], "CHUCK_NORRIS_EXCEPTION: Chuck Norris can dereference NULL.", max_quote_size);

  strncpy(quote_list[29],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s code doesn’t follow a coding convention. It is the coding convention.",
        max_quote_size);
  strncpy(quote_list[30],
        "CHUCK_NORRIS_EXCEPTION: Users don’t mark Chuck Norris’s answers as accepted. The universe accepts them out of a sense of truth and justice.",
        max_quote_size);
  strncpy(quote_list[31],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris can divide by Zero.",
        max_quote_size);
  strncpy(quote_list[32],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s Stack Overflow reputation is only as modest as it is because of integer overflow (SQL Server does not have a datatype large enough).",
        max_quote_size);
  strncpy(quote_list[33],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris is the only top 100 Stack Overflow user who is human. The others are bots that he coded to pass the time between questions.",
        max_quote_size);
  strncpy(quote_list[34],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris coded his last project entirely in Microsoft Paint, just for the challenge.",
        max_quote_size);
  strncpy(quote_list[35],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris does not use exceptions when programming. He has not been able to identify any of his code that is not exceptional.",
        max_quote_size);
  strncpy(quote_list[36],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris’s code fails to compile, the compiler apologizes.",
        max_quote_size);
  strncpy(quote_list[37],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris does not use revision control software. None of his code has ever needed revision.",
        max_quote_size);
  strncpy(quote_list[38],
        "CHUCK_NORRIS_EXCEPTION: When you search for “guru” on Google it says “Did you mean Chuck Norris?”",
        max_quote_size);
  strncpy(quote_list[39],
        "CHUCK_NORRIS_EXCEPTION: There are two types of programmers: good programmers, and those who are not Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[40],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris points to null, null quakes in fear.",
        max_quote_size);
  strncpy(quote_list[41],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris is the traveling salesman. Only he knows the shortest route.",
        max_quote_size);
  strncpy(quote_list[42],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris took the red pill and the blue pill, and can phase-shift in and out of the Matrix at will.",
        max_quote_size);
  strncpy(quote_list[43],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris pushes a value onto a stack, it stays pushed.",
        max_quote_size);
  strncpy(quote_list[44],
        "CHUCK_NORRIS_EXCEPTION: When invoking one of Chuck Norris’s callbacks, the runtime adds “please”.",
        max_quote_size);
  strncpy(quote_list[45],
        "CHUCK_NORRIS_EXCEPTION: Drivers think twice before they dare interrupt Chuck Norris’s code.",
        max_quote_size);
  strncpy(quote_list[46],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris does not sleep… He waits.",
        max_quote_size);
  strncpy(quote_list[47],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t answer questions on Stack Overflow... He stares them down till they answer themselves.",
        max_quote_size);
  strncpy(quote_list[48],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris can stop an infinite loop just by thinking about it.",
        max_quote_size);
  strncpy(quote_list[49],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t need a debugger, he just stares down the bug until the code confesses.",
        max_quote_size);
  strncpy(quote_list[50],
        "CHUCK_NORRIS_EXCEPTION: There is no ‘CTRL’ button on Chuck Norris’s computer. Chuck Norris is always in control.",
        max_quote_size);
  strncpy(quote_list[51],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris won the “Hello World in less than 20 bytes” contest by developing a single byte program.",
        max_quote_size);
  strncpy(quote_list[52],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris does not resolve software problems. The problems resolve themselves the moment he walks into the office.",
        max_quote_size);
  strncpy(quote_list[53],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris can answer a question well before it is asked and then get several up-votes whilst he has yet to finish typing the solution.",
        max_quote_size);
  strncpy(quote_list[54],
        "CHUCK_NORRIS_EXCEPTION: The Chuck Norris badge is awarded for posting a better answer than Chuck Norris. Only Chuck Norris can earn this badge.",
        max_quote_size);
  strncpy(quote_list[55],
        "CHUCK_NORRIS_EXCEPTION: God said: ‘Let there be light,’ only so he could see what Chuck Norris was up to.",
        max_quote_size);
  strncpy(quote_list[56],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s keyboard doesn’t have F1 key, the computer asks for help from him.",
        max_quote_size);
  strncpy(quote_list[57],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris presses Ctrl+Alt+Delete, a worldwide computer restart is initiated. The same goes for format.",
        max_quote_size);
  strncpy(quote_list[58],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris uses Visual Studio to burn CDs.",
        max_quote_size);
  strncpy(quote_list[59],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris is not close to perfection, perfection is close to Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[60],
        "CHUCK_NORRIS_EXCEPTION: God didn’t really create the world in six days, because Chuck Norris optimized it to one.",
        max_quote_size);
  strncpy(quote_list[61],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris dreams in ones and zeros. When two shows up, it is a nightmare. But again that’s only in theory. Two doesn’t exist for Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[62],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris solves an equation, the variables become constants.",
        max_quote_size);
  strncpy(quote_list[63],
        "CHUCK_NORRIS_EXCEPTION: If anyone writes delete ChuckNorris; in C, the Apocalypse will come.",
        max_quote_size);
  strncpy(quote_list[64],
        "CHUCK_NORRIS_EXCEPTION: Once Chuck Norris went to the library… Since then the library was dynamically linked.",
        max_quote_size);
  strncpy(quote_list[65],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris has the key to Open Source. He just doesn’t want to close it.",
        max_quote_size);
  strncpy(quote_list[66],
        "CHUCK_NORRIS_EXCEPTION: Compatibility doesn’t exist in Chuck Norris’s dictionary. He can easily work in Microsoft Office in Linux on a Mac.",
        max_quote_size);
  strncpy(quote_list[67],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris is programming, the Garbage Collector rests. The objects know when to destroy themselves.",
        max_quote_size);
  strncpy(quote_list[68],
        "CHUCK_NORRIS_EXCEPTION: If the Internet is the web, then Chuck Norris is the spider.",
        max_quote_size);
  strncpy(quote_list[69],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris is on a diet and doesn’t eat fast food, all hard disks change from FAT to NTFS.",
        max_quote_size);
  strncpy(quote_list[70],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris has written the best programming language. Its source has just one command… void Chuck NorrisSkeet();",
        max_quote_size);
  strncpy(quote_list[71],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t use #include. He thinks of it as cheating.",
        max_quote_size);
  strncpy(quote_list[72],
        "CHUCK_NORRIS_EXCEPTION: When a null reference exception goes to sleep, it checks under the bed for Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[73],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t need delegates, he does all the work himself.",
        max_quote_size);
  strncpy(quote_list[74],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t call a background worker, background workers call Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[75],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t write books, the words assemble themselves out of fear.",
        max_quote_size);
  strncpy(quote_list[76],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris throws an exception, nothing can catch it.",
        max_quote_size);
  strncpy(quote_list[77],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris is beyond Turing-complete; he is Turing-invincible.",
        max_quote_size);
  strncpy(quote_list[78],
        "CHUCK_NORRIS_EXCEPTION: There is simply no Halting Problem within a 10-meter radius of John Skeet, because computers ALWAYS halt in his presence.",
        max_quote_size);
  strncpy(quote_list[79],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t look for reputation. Reputation looks for Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[80],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris can do pair programming with himself.",
        max_quote_size);
  strncpy(quote_list[81],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris installed Visual Studio he opted not to install the debugger.",
        max_quote_size);
  strncpy(quote_list[82],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris saves a file, the file thanks him.",
        max_quote_size);
  strncpy(quote_list[83],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris is immutable. If something’s going to change, it’s going to have to be the rest of the universe.",
        max_quote_size);
  strncpy(quote_list[84],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s addition operator doesn’t commute; it teleports to where he needs it to be.",
        max_quote_size);
  strncpy(quote_list[85],
        "CHUCK_NORRIS_EXCEPTION: Anonymous methods and anonymous types are really all called Chuck Norris. They just don’t like to boast.",
        max_quote_size);
  strncpy(quote_list[86],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris doesn’t have performance bottlenecks. He just makes the universe wait its turn.",
        max_quote_size);
  strncpy(quote_list[87],
        "CHUCK_NORRIS_EXCEPTION: Jeff Atwood bought a monster GPU just to calculate Chuck Norris’s rep on Stack Overflow. CPUs don’t cut it anymore.",
        max_quote_size);
  strncpy(quote_list[88],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris does a search on Google... the only result is “I’ll be right back”.",
        max_quote_size);
  strncpy(quote_list[89],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris returned IntelliSense and got his money back!",
        max_quote_size);
  strncpy(quote_list[90],
        "CHUCK_NORRIS_EXCEPTION: When Chuck Norris presses F5, the Garbage collector collects itself.. there is no other garbage.",
        max_quote_size);
  strncpy(quote_list[91],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris once wrote an entire operating system in his sleep on a Treo with no battery, powered only by the force of his will.",
        max_quote_size);
  strncpy(quote_list[92],
        "CHUCK_NORRIS_EXCEPTION: The only time Chuck Norris was wrong was when he thought he had made a mistake.",
        max_quote_size);
  strncpy(quote_list[93],
        "CHUCK_NORRIS_EXCEPTION: If you have 10000 reputation points and Chuck Norris has 10000 reputation points, Chuck Norris has more reputation than you.",
        max_quote_size);
  strncpy(quote_list[94],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris does not run his programs. He just whispers “you better run”. And it runs.",
        max_quote_size);
  strncpy(quote_list[95],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris was once second in rank, behind Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[96],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris codes only with final sealed methods. No one has ever needed to override any of Chuck Norris’s code.",
        max_quote_size);
  strncpy(quote_list[97],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris is IntelliSense.",
        max_quote_size);
  strncpy(quote_list[98],
        "CHUCK_NORRIS_EXCEPTION: Chuck Norris’s heart rate is 5 GHz.",
        max_quote_size);
  strncpy(quote_list[99],
        "CHUCK_NORRIS_EXCEPTION: Private methods in other libraries become public automatically once required in Chuck Norris’s code.",
        max_quote_size);
  strncpy(quote_list[100],
        "CHUCK_NORRIS_EXCEPTION: When Yoda needs advice, he calls Chuck Norris.",
        max_quote_size);
  strncpy(quote_list[101],
        "CHUCK_NORRIS_EXCEPTION: Only Chuck Norris earned the coveted “Chuck Norris” badge.",
        max_quote_size);
  strncpy(quote_list[102],
        "CHUCK_NORRIS_EXCEPTION: Norrised: The act of attempting to answer a Stack Overflow question only to find out that Chuck Norris has already answered it definitively and better than you would have ever done.",
        max_quote_size);
  strncpy(quote_list[103],
        "CHUCK_NORRIS_EXCEPTION: If Chuck Norris posts a duplicate question on StackOverflow, the original question will be closed as a duplicate.",
        max_quote_size);
        
  ind = (int)((double)rand() * (104.0 + 1.0) / RAND_MAX);

  fprintf(stderr, "%s\n", quote_list[ind]);
  fflush(stderr);
}

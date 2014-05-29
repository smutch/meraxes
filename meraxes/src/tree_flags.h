#define TREE_CASE_NO_PROGENITORS                1         // Set for halos that have no progenitors.
#define TREE_CASE_MAIN_PROGENITOR               2         // Set for the progenitor with the highest match score. (propagated for ghosts)
#define TREE_CASE_STRAYED                       4         // Set for halos for which a descendant was not found
#define TREE_CASE_REMNANT                       8         // Set for halos with more than one progenitor.
#define TREE_CASE_MERGER                        16        // Set when new IDs are created (ie. last point the halo was seen).
                                                          //    Set only for the last ghost in ghost-populated trees for mergers w/ offset>1.
#define TREE_CASE_DROPPED                        32       // Set if file_offset>1 and TREE_CASE_MATCHED_TO_BRIDGE is not set
#define TREE_CASE_BRIDGED                        64       // Set for halos with multiple unique back-matches from halos with unique IDs
#define TREE_CASE_EMERGED                        128      // Set when a match is made identifying this halo as emerged
#define TREE_CASE_FRAGMENTED_NEW                 256      // Set for halos that have been marked TREE_CASE_EMERGED_CANDIDATE but not TREE_CASE_EMERGED
                                                          //    (unless it's the backmatch with the most massive descendant; that halo is considered
                                                          //     to be the source of any fragmented halos)
#define TREE_CASE_FRAGMENTED_STRAYED              512     // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                          //    decendant_id!=a valid id (ie they are not a progenitor of anything). (propagated for ghosts)
#define TREE_CASE_FRAGMENTED_RETURNED             1024    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                          //    decendant_id==the id of the halo they are emerged from. (propagated for ghosts)
#define TREE_CASE_FRAGMENTED_EXCHANGED            2048    // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose 
                                                          //    decendant_id!=the id of the halo they are emerged but is nevertheless valid 
                                                          //    (ie. they are still a progenitor of something). (propagated for ghosts)
#define TREE_CASE_EMERGED_CANDIDATE               4096    // Set when a halo is identified as a unique back-match to a halo marked TREE_CASE_BRIDGED 
                                                          //    and is not identified as the BRIDGE's main descendant
#define TREE_CASE_MATCHED_TO_EMERGED              8192    // Set when a halo is matched to an emerged halo
#define TREE_CASE_GHOST                           16384   // Marks ghost halos in ghost-populated trees
#define TREE_CASE_GHOST_NULL                      32768   // Marks a ghost halo where a subgroup is it's own group.
                                                          //    This is a default behaviour that occurs when a group is strayed but one of 
                                                          //    it's subgroups isn't.
#define TREE_CASE_UNPROCESSED                     65536   // For internal use.  This should never be seen in the output.
#define TREE_CASE_INVALID                         131072  // For internal use.  This should never be seen in the output.


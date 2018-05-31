// TTPXX=Two-to-the-power-XX
#define TTTP00 1
#define TTTP01 2
#define TTTP02 4
#define TTTP03 8
#define TTTP04 16
#define TTTP05 32
#define TTTP06 64
#define TTTP07 128
#define TTTP08 256
#define TTTP09 512
#define TTTP10 1024
#define TTTP11 2048
#define TTTP12 4096
#define TTTP13 8192
#define TTTP14 16384
#define TTTP15 32768
#define TTTP16 65536
#define TTTP17 131072
#define TTTP18 262144
#define TTTP19 524288
#define TTTP20 1048576
#define TTTP21 2097152
#define TTTP22 4194304
#define TTTP23 8388608
#define TTTP24 16777216
#define TTTP25 33554432
#define TTTP26 67108864
#define TTTP27 134217728
#define TTTP28 268435456
#define TTTP29 536870912
#define TTTP30 1073741824
#define TTTP31 2147483648
#define TTTP32 4294967296

#define TREE_CASE_NO_PROGENITORS TTTP00 // Set for halos that have no progenitors.
#define TREE_CASE_MAIN_PROGENITOR TTTP01 // Set for the progenitor with the highest match score.
#define TREE_CASE_MOST_MASSIVE TTTP02 // Marks the most massive substructure.
#define TREE_CASE_DOMINANT TTTP03 // Marks the dominant     substructure.
#define TREE_CASE_REMNANT TTTP04 // Set for halos with more than one progenitor.
#define TREE_CASE_MERGER_PRIMARY TTTP05 // Set when a halo is deemed to be the primary   progenitor of a merger
#define TREE_CASE_MERGER TTTP06 // Set when a halo is deemed to be the secondary progenitor of a merger
#define TREE_CASE_STRAYED TTTP07 // Set for halos for which a descendant was not found
#define TREE_CASE_DROPPED TTTP08 // Set if file_offset>1 and TREE_CASE_MATCHED_TO_BRIDGE is not set
#define TREE_CASE_BRIDGED TTTP09 // Set for halos with multiple unique back-matches from halos with unique IDs
#define TREE_CASE_EMERGED TTTP10 // Set when a match is made identifying this halo as emerged
#define TREE_CASE_FRAGMENTED_NEW TTTP11 // Set for halos that have been marked TREE_CASE_EMERGED_CANDIDATE but not TREE_CASE_EMERGED
//    (unless it's the backmatch with the most massive descendant; that halo is considered
//     to be the source of any fragmented halos)
#define TREE_CASE_FRAGMENTED_STRAYED TTTP12 // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), and whose
//    decendant_id!=a valid id (ie they are not a progenitor of anything).
#define TREE_CASE_FRAGMENTED_NORMAL TTTP13 // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), but
//    are otherwise part of a connected progenitor line
#define TREE_CASE_FRAGMENTED_EJECTED TTTP14 // Set for halos that are marked TREE_CASE_FRAGMENTED (see below), for which
//    a progenitor was possible but could not be used because it was needed
//    as a progenitor of the bridged halo this fragmented halo is back matched to.
#define TREE_CASE_EMERGED_CANDIDATE TTTP15 // Set when a halo is identified as a unique back-match to a halo marked TREE_CASE_BRIDGED and
//    is not selected as the bridged halo's descendant.  This switch is turned off if the halo
//    is marked TREE_CASE_EMERGED or TREE_CASE_FRAGMENTED_NEW.
#define TREE_CASE_MATCHED_TO_EMERGED TTTP16 // Set when a halo is matched to an emerged halo
#define TREE_CASE_2WAY_MATCH TTTP17 // Set when the match between a halo and it's descendant is mutual
#define TREE_CASE_SET_BY_BACKMATCH TTTP18 // Set in cases where a descendant was set using backmatch information, rather
//    than forematch information.  In these cases, the halo's given matching score and
//    2way flag are from the back match info, not the forematch info.
#define TREE_CASE_DOMINANT_DROPPED TTTP19 // Marks a group whose dominant substructure is presently dropped
#define TREE_CASE_GHOST TTTP20 // Marks ghost halos in ghost-populated trees
#define TREE_CASE_GHOST_NULL TTTP21 // Marks a ghost halo where a subgroup is it's own group.
//    This is a default behaviour that occurs when a group is strayed but one of
//    it's subgroups isn't.
#define TREE_CASE_HIGHRES TTTP22 // Marks a halo connected from the highres simulation.
#define TREE_CASE_HIGHRES_TRAITOR TTTP23 // Marks a halo connected from the highres simulation which does not have any descendants
//    because it merges into another group in the highres
#define TREE_CASE_HIGHRES_CONNECTION TTTP24 // Marks a halo connected from the highres simulation to a halo in the lowres sim
#define TREE_CASE_HIGHRES_HERMIT TTTP25 // Marks a halo connected from the highres simulation whose descendant cannot be resolved in lowres

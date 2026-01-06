#define _GNU_SOURCE
#include "dsm.h"
#include <stdio.h>
#include <stdlib.h>

// Mock sa_strmap if needed, or link against sa.o
// sa_strmap is in sa.c. We can link against it.

char *base_chars = "-AGCUN"; // 0=Gap, 1=A, ...

int main() {
  printf("Dumping DSM_T04_POS [q1][q2][t1][t2] (using GAP=%d)\n", GAP);

  // Iterate 0..5 for indices
  for (int q1 = 0; q1 < 6; q1++) {
    for (int q2 = 0; q2 < 6; q2++) {
      for (int t1 = 0; t1 < 6; t1++) {
        for (int t2 = 0; t2 < 6; t2++) {
          // dsm_t04_pos is defined in dsm.h/dsm-t04.c
          int val = dsm_t04_pos[q1][q2][t1][t2];
          if (val != -2000 && val != 0 && val != -123 && val != -22) {
            // Print non-trivial values to check scale
          }
          // Print header for first few
          if (q1 == 0 && q2 == 1 && t1 == 0 && t2 == 0) { // Gap-A vs Gap-Gap
            printf("[Gap][A][Gap][Gap] = %d\n", val);
          }
        }
      }
    }
  }

  // Specific check for the problematic case?
  // Gap Open costs
  // Q: Gap-A (Gap at 0, A at 1)
  // T: Gap-Gap

  // Let's just print the whole table in a format Rust can paste?
  // Or just print indices to confirm mapping.

  printf("pub const DSM_T04_POS_DUMP: DsmTable = [\n");
  for (int i = 0; i < 6; i++) {
    printf("    [\n");
    for (int j = 0; j < 6; j++) {
      printf("        [\n");
      for (int k = 0; k < 6; k++) {
        printf("            [");
        for (int l = 0; l < 6; l++) {
          printf("%d, ", dsm_t04_pos[i][j][k][l]);
        }
        printf("],\n");
      }
      printf("        ],\n");
    }
    printf("    ],\n");
  }
  printf("];\n");

  return 0;
}

/**
 * output_utils.c
 *
 * Shared utility functions for HDF5 output operations.
 */

#include "output_utils.h"
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>

int
read_env_int_or_default(const char *name, int def) {
  const char *v = getenv(name);
  return v ? atoi(v) : def;
}

int
ensure_directory_exists(const char *dir) {
  char   path[PATH_MAX];
  size_t len = strlen(dir);

  if (len == 0 || len >= sizeof(path))
    return -1;

  snprintf(path, sizeof(path), "%s", dir);

  if (path[len - 1] == '/' && len > 1)
    path[len - 1] = '\0';

  for (char *p = path + 1; *p != '\0'; ++p) {
    if (*p != '/')
      continue;
    *p = '\0';
    if (mkdir(path, 0777) < 0 && errno != EEXIST)
      return -1;
    *p = '/';
  }

  if (mkdir(path, 0777) < 0 && errno != EEXIST)
    return -1;

  return 0;
}

void
build_output_filepath(char *buf, size_t sz, const char *file_name) {
  const char *base_dir_env = getenv("HDF_OUTPUT_DIR");

  const char *base_dir = NULL;

  if (base_dir_env == NULL) {
    base_dir = "./";
  } else {
    base_dir = base_dir_env;
  }

  char dir_name[] = "3d_field";

  char full_path[PATH_MAX];
  snprintf(full_path, sizeof(full_path), "%s/%s", base_dir, dir_name);

  if (ensure_directory_exists(full_path) < 0) {
    fprintf(stderr, "[output] failed to create directory: %s\n", full_path);
    if (sz > 0)
      buf[0] = '\0';
    return;
  }

  if (file_name == NULL || *file_name == '\0') {
    snprintf(buf, sz, "%s", full_path);
    return;
  }

  snprintf(buf, sz, "%s/%s", full_path, file_name);
}

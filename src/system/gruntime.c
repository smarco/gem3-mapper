/*
 *  GEM-Mapper v3 (GEM3)
 *  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
 *
 *  This file is part of GEM-Mapper v3 (GEM3).
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * PROJECT: GEM-Mapper v3 (GEM3)
 * AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
 */

#include "system/gruntime.h"
#include "system/errors.h"
#include "system/mm.h"
#include "system/gthread.h"
#include "profiler/profiler_gem.h"

/*
 * GEM Runtime
 */
void gruntime_init(
    const uint64_t num_threads,
    char* const tmp_folder) {
  // GEM error handler
  gem_handle_error_signals();
  // Setup Profiling (Add master thread)
  PROF_NEW(num_threads);
  // Register Master-Thread
  gem_thread_register_id(0);
  // Setup temporal folder
  if (tmp_folder!=NULL) mm_set_tmp_folder(tmp_folder);
}
void gruntime_destroy(void) {
  // Clean-up Profiler
  PROF_DELETE();
  // Clean-up Thread Info
  gem_thread_cleanup();
}


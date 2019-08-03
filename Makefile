###############################################################################
#  GEM-Mapper v3 (GEM3)
#  Copyright (c) 2011-2017 by Santiago Marco-Sola  <santiagomsola@gmail.com>
#
#  This file is part of GEM-Mapper v3 (GEM3).
#
#  This program is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# PROJECT: GEM-Mapper v3 (GEM3)
# AUTHOR(S): Santiago Marco-Sola <santiagomsola@gmail.com>
###############################################################################

# Definitions
ROOT_PATH=$(CURDIR)
include Makefile.mk

all: release

# Optimized 
release: setup
	$(MAKE) --directory=resources release
	$(MAKE) --directory=src release
	
# Optimized + DebugSymbols
devel: setup
	$(MAKE) --directory=resources devel
	$(MAKE) --directory=src devel
	
# Optimized + DebugSymbols + GEMProfile
profile: setup
	$(MAKE) --directory=resources profile
	$(MAKE) --directory=src profile
	
# DebugSymbols + GEMProfile + GEMDebug
debug: setup
	$(MAKE) --directory=resources debug
	$(MAKE) --directory=src debug
	
complexity:
	pmccabe -vt include/*/*.h src/*/*.c
	
setup: 
	@-mkdir -p $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	@-ln -s $(FOLDER_RESOURCES) $(FOLDER_INCLUDE)/resources 2> /dev/null ||:

clean:
	$(MAKE) --directory=resources clean
	@rm -rf $(FOLDER_BIN) $(FOLDER_BUILD) $(FOLDER_LIB)
	@rm -rf $(FOLDER_INCLUDE)/resources

distclean: clean
	@rm -f config.status config.log Makefile.mk

archive:
	# Taken from https://ttboj.wordpress.com/2015/07/23/git-archive-with-submodules-and-tar-magic/
	@echo Running git archive...
	# Use HEAD if tag doesn't exist yet, so that development is easier...
	git archive --prefix=gem3-mapper/ -o $(ARCHIVENAME) $(VERSION) 2> /dev/null || (echo 'Warning: $(VERSION) does not exist.' && git archive --prefix=gem3-mapper/ -o $(ARCHIVENAME) HEAD)
	# TODO: if git archive had a --submodules flag this would easier!
	mkdir -p archivebuild
	@echo Running git archive submodules...
	# I thought I would need --ignore-zeros, but it doesn't seem necessary!
	p=`pwd` && (echo .; git submodule foreach) | while read entering path; do \
	    temp="$${path%\'}"; \
	    temp="$${temp#\'}"; \
	    path=$$temp; \
	    [ "$$path" = "" ] && continue; \
	    (cd $$path && git archive --prefix=gem3-mapper/$$path/ HEAD > $$p/archivebuild/tmp.tar && tar --concatenate --file=$$p/$(ARCHIVENAME) $$p/archivebuild/tmp.tar && rm $$p/archivebuild/tmp.tar); \
	done
	rm -rf archivebuild
	gzip $(ARCHIVENAME)


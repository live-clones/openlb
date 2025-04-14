###########################################################################
## Makefile for OpenLB's automatic CSE workflow

# Directory containing the target files with dynamics type information
TARGET_DIR = targets/dynamics

# List of targets to extract expressions
EXT_TARGET_LIST = target.txt

# Define the directory containing the template for extractions
EXT_DIR = extractions

# Define the directory containing the .out files
EXP_DIR = expressions

# Define the directory for dumping the .cse.h files
CSE_DIR = include

# File for including the generated .cse.h files
INCLUDE_HEADER := $(EXP_DIR)/$(CSE_DIR)/generated_cse.h

##########################################################################
.PHONY: target_list
target_list:
	bash unite_targets.sh $(TARGET_DIR) target_tmp.txt;
	bash format_targets.sh target_tmp.txt $(EXT_TARGET_LIST);
	rm target_tmp.txt;

##########################################################################
.PHONY: decompose_targets
decompose_targets: $(target_list)
	mkdir -p $(EXP_DIR)/$(CSE_DIR)/; \
	bash decompose_targets.sh $(EXT_TARGET_LIST) $(EXT_DIR) $(INCLUDE_HEADER)

EXT_TARGETS := $(wildcard $(EXT_DIR)/*.txt)

EXT_OBJECTS := $(EXT_TARGETS:.txt=.out)

$(EXT_DIR)/%.out: $(EXT_DIR)/%.txt
	@mkdir -p $(EXT_DIR)/$*;
	@if [ ! -f "$(EXP_DIR)/$*.out" ]; then \
		content=$$(cat $<); \
		cd $(EXT_DIR)/$*; \
		cp ../extract.cpp.template ../Makefile .; \
		python ../../extraction.py "$$content" "extract.cpp.template" "$<" "./" "$(EXP_DIR)/" ; \
	else \
		echo "Skip existing $(EXP_DIR)/$*.out."; \
	fi

.PHONY: check_targets
check_targets: $(EXT_OBJECTS)
	bash check_extraction.sh

##########################################################################
EXP_SOURCES := $(wildcard $(EXP_DIR)/*.out)

CSE_OBJECTS := $(EXP_SOURCES:.out=.cse.h)

$(EXP_DIR)/%.cse.h: $(EXP_DIR)/%.out
	@mkdir -p $(EXP_DIR)/$(CSE_DIR);
	@if [ ! -f "$(EXP_DIR)/$(CSE_DIR)/$*.cse.h" ]; then \
		python optimize_cse.py "$<" "$(EXP_DIR)/dynamics.cse.template" "$(EXP_DIR)/$(CSE_DIR)/$*.cse.h"; \
	else \
		echo "Skip existing $(EXP_DIR)/$(CSE_DIR)/$*.cse.h."; \
	fi

.PHONY: check_optimized
check_optimized: $(CSE_OBJECTS)
	bash check_optimized.sh
##########################################################################

prepare_extraction: target_list decompose_targets

extraction: $(EXT_OBJECTS) check_targets

optimize_cse: $(CSE_OBJECTS) check_optimized

##########################################################################

.PHONY: clean_cse
clean_cse:
	rm -f $(EXP_DIR)/$(CSE_DIR)/*.cse.h $(INCLUDE_HEADER)

.PHONY: clean_extraction
clean_extraction:
	rm -f $(EXP_DIR)/*.out $(EXT_DIR)/*.txt
	@for dir in $(EXT_DIR)/*; do\
		if [ -d "$$dir" ]; then \
			rm -rf "$$dir"/; \
		fi \
	done

.PHONY: clean_target_list
clean_target_list:
	rm $(EXT_TARGET_LIST)

.PHONY: clean
clean: clean_extraction

.PHONY: clean_all
clean_all: clean_cse clean_extraction clean_target_list


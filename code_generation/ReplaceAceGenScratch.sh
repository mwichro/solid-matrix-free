
#!/bin/bash

# Check if file is provided as argument
if [ -z "$1" ]; then
  echo "Usage: $0 <source_file>"
  exit 1
fi

echo -e "\e[31mHello, World!\e[0m"

SOURCE_FILE=$1
TARGET_FILE="${SOURCE_FILE%.cpp}.cc"
echo ---------


echo -e "Reading \e[31m $SOURCE_FILE \e[0m"
echo -e "Saving modified  to   \e[32m $TARGET_FILE \e[0m"

TEMP_VARS=()


cp $SOURCE_FILE $TARGET_FILE


# Function to run clang-format on a C++ source file with a specific .clang-format file
run_clang_format() {


  # Get the directory of the script
  script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

  # Path to the .clang-format file
  clang_format_file="$script_dir/../../.clang-format"

  local file=$1

  clang-format -i -style=file "$file" --assume-filename="$clang_format_file"
  #echo "Ran clang-format on $file with $clang_format_file"
}



# Function to remove double acegen_scratch[] from the function header
remove_acegen_scratch() {
   
  local file=$1
  local first_brace_line=$(sed -n '/{/=' "$file" | head -n 1)
  echo Removing acegen_scratch from the header
  echo Function body starts in line: $first_brace_line
  sed -i "1,${first_brace_line}  s/double *acegen_scratch\[[0-9]*\],* *//g" "$file"
}


# Function to replace acegen_scratch[XXXX] with acegen_scratch__XXXX__
replace_and_collect() {
  local file=$1
  # Use sed to find and replace while collecting unique temporary variables
  local replaced_content=$(sed -E 's/acegen_scratch\[([0-9]+)\]/acegen_scratch__\1__/g' "$file")
  
  # Extract unique temporary variables
  local vars=$(echo "$replaced_content" | grep -oP 'acegen_scratch__[0-9]+__' | sort | uniq)
  
  # Save replaced content back to the file
  echo "$replaced_content" > "$file"
  
  # Collect temporary variables
  for var in $vars; do
    TEMP_VARS+=("$var")
  done

  # Get the line number of the first { and quit processing
  local first_brace_line=$(sed -n '/{/=' "$file" | head -n 1)

  # Generate the line to declare all temporary variables
  declare_line="Number $(IFS=,; echo "${TEMP_VARS[*]}");"

  # Add empty lines and a comment before and after the added content
  content_to_add="
$declare_line
"
  #echo ---------
  #echo $content_to_add
  #echo ---------


  # Add the content after the first { in the file
  sed -i "${first_brace_line}a\
  \\$content_to_add
  " "$file"


  echo "--Added line to declare all temporary variables"

}

# Function to add template<typename Number> before the first occurrence of void in the file
add_template_before_void() {
  local file="$1"

  # Get the line number of the first occurrence of void
  local first_void_line=$(grep -n '^void' "$file" | head -n 1 | cut -d ':' -f 1)

  # Add template<typename Number> before the first occurrence of void
  sed -i "${first_void_line}i\\
template<typename Number>
" "$file"

  echo "Added template<typename Number> before the first occurrence of void in $file"
}



# Function to replace double [] and double [][] with Tensor types
replace_double_with_tensor() {
  local file=$1
  
  local first_brace_line=$(sed -n '/{/=' "$file" | head -n 1)
  # Replace double [][] with Tensor<2, dim, Number> within the first few lines
  sed -i "1,${first_brace_line} s/double \([a-zA-Z_][a-zA-Z0-9_]*\)\[[0-9]\+\]\[[0-9]\+\]/Tensor<2, dim, Number> \& \1/g" "$file"


  # Replace double [] with Tensor<1, dim, Number> within the first few lines
  sed -i "1,${first_brace_line} s/double \([a-zA-Z_][a-zA-Z0-9_]*\)\[[0-9]\+\]/Tensor<1, dim, Number> \& \1/g" "$file"

}


# Run clang-format on the specified source file
run_clang_format "$TARGET_FILE"

# Remove double acegen_scratch[] from the header
remove_acegen_scratch "$TARGET_FILE"
run_clang_format "$TARGET_FILE"

# Replace and collect temporary variables
replace_and_collect "$TARGET_FILE"
run_clang_format "$TARGET_FILE"

#add_template_before_void "$TARGET_FILE"

# Replace double types with Tensor types in the source file
#replace_double_with_tensor "$TARGET_FILE"
#run_clang_format "$TARGET_FILE"


# Generate the line to declare all temporary variables
#declare_line="Number $(IFS=,; echo "${TEMP_VARS[*]}");"

# Output the line to declare all temporary variables
#echo "$declare_line"
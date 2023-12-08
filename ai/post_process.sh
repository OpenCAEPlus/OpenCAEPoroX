#!/bin/bash

# Define source and destination folders
source_folder="doc_include"
bad_folder="bad"
good_folder="good"
manual_folder="manual"
start_pattern="^\`\`\`"
end_pattern="^\`\`\`$"

# Check if source folder exists
if [[ ! -d "$source_folder" ]]; then
  echo "Error: Source folder '$source_folder' does not exist."
  exit 1
fi

# Check if destination folder exists
create_folder() {
  local folder_name="$1"
  if [[ ! -d "$folder_name" ]]; then
    echo "Creating destination folder '$folder_name'..."
    mkdir -p "$folder_name" || exit 1
  fi
}

create_folder "$bad_folder"
create_folder "$good_folder"
create_folder "$manual_folder"

# Loop through all files in the source folder
for file in "$source_folder"/*; do
  # Check if file is empty
  if [[ ! -s "$file" ]]; then
    echo "Moving empty file '$file' to '$bad_folder'..."
    cp "$file" "$bad_folder" || exit 1
  else
    # echo "nothing"
    # Check for the presence of patterns
    end_count=$(grep -c "$end_pattern" "$file")

    # Check if the ending pattern is present
    if [[ $end_count -ne 1 ]]; then
      echo "Moving incomplete file '$file' to '$manual_folder'..."
      cp "$file" "$manual_folder" || exit 1
    else
      # Get the content between ``` lines

      # Remove any lines before or after the content
      # content=$(sed -n "2,$(wc -l <<< \"$content\")p" <<< "$content")
      # Determine target folder
      cp "$file" "$good_folder"/ || exit 1
      echo "Moved '$file' to '$good_folder'."
      sed -i -n -e "/$start_pattern/,/$end_pattern/p" $good_folder/$(basename "$file")
      sed -i '1d;$d' $good_folder/$(basename "$file")
    fi
  fi
done

echo "Done! Processed all files in '$source_folder'."

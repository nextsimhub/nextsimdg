import os
import re
import string

# Prerequisites:
# 
# github_changelog_generator repo -t token --date-format '%Y-%m-%d %H:%M' -o 'CHANGELOG.md'
# git tag -l -n9 > a.txt
# All tags should be of type 'v\d[\.\w]*'

if __name__ == '__main__':
    with open('a.txt', 'r') as foo:
        data = foo.readlines()
        foo.close()
    tags_title_annotation = {}
    matched = None
    for line in data:
        if re.search('v\d[\.\w]*', line):
            matched = re.findall('v\d[\.\w]*', line)[0]
            tags_title_annotation[matched] = {
                'title': '',
                'annotation': []
            }
            title = line.replace(matched, '').strip()
            title = title[0].upper() + title[1:]
            tags_title_annotation[matched]['title'] = title
        elif matched is not None:
            line = line.strip()
            annotation = line.replace(matched, '')
            tags_title_annotation[matched]['annotation'].append(annotation)

    with open('CHANGELOG.md', 'r') as fil:
        data = fil.readlines()
        fil.close()
    matched = None
    data_w = data
    for index, line in enumerate(data):
        if re.search('## \[v\d[\.\w]*\]', line):
            matched = re.findall('v\d[\.\w]*', line)[0]

            insert_data = [
                "**Description:**",
                tags_title_annotation[matched]['title']
            ]
            values = "\n" + \
                "\n".join(insert_data) + \
                "\n" + \
                "\n".join(tags_title_annotation[matched]['annotation']) + \
                "\n"
            data_w.insert(index + 1, values)
    
    with open('CHANGELOG.md', 'w') as fin:
        fin.writelines(data_w)
        fin.close()


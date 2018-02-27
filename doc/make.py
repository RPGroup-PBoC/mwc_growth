#! /usr/bin/env python
import os
import frontmatter


def main():
    with open('build.yaml') as f:
        m, _ = frontmatter.parse(f.read())

        code = """
        pandoc default.yaml -i {} --filter=pandoc-eqnos --filter=pandoc-crossref -o {}.pdf
        """.format(m['include'], m['name'])
        os.system(code)


if __name__ == '__main__':
    main()
    print('document successfully compiled')

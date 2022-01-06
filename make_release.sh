#!/bin/sh
set -e

export PARSED_VERSION=$(head -n1 src/pygcgopt/__init__.py | cut -d "'" -f2)

echo "Parsed version number ${PARSED_VERSION} from src/pygcgopt/__init__.py."

if ! grep -q "v${PARSED_VERSION}" CHANGELOG.md; then
    echo "The parsed version number ${PARSED_VERSION} is not contained in CHANGELOG.md. Please update the changelog and try again."
    exit 1
fi

echo "Updating remotes"
git remote update

if [[ `git status --porcelain` ]]; then
    echo "Your git working tree is not clean. Please commit all changes and push them before crafting the release."
    exit 1
fi

git push

echo "Creating release tag"
gh release create "v${PARSED_VERSION}" --title "v${PARSED_VERSION}" --notes "Release ${PARSED_VERSION}"


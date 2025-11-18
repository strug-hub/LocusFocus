#!/bin/bash
# Sets up the locusfocus user and group. Run with sudo. needed for docker permissions

set -euox pipefail

if [[ -z $1 ]]; then
    echo "No username provided. Please provide your current username."
    echo "Usage: setup-locusfocus-user.sh <username>"
    exit 1
fi

# Add group if it doesn't exist
if [[ -z $(getent group locusfocus) ]]; then
    groupadd --gid 4644 locusfocus
fi

# Add user if it doesn't exist
if [[ -z $(getent passwd locusfocus) ]]; then
    useradd locusfocus -u 4644 -g 4644 -m -s /bin/bash
fi

# Update file permissions

chown -R locusfocus:locusfocus ./app
chown -R locusfocus:locusfocus ./tests
chown -R locusfocus:locusfocus ./README.md
chown -R locusfocus:locusfocus ./pyproject.toml
chown -R locusfocus:locusfocus ./poetry.lock
chown -R locusfocus:locusfocus ./data/cache

chmod -R g=u ./app
chmod -R g=u ./tests
chmod -R g=u ./README.md
chmod -R g=u ./pyproject.toml
chmod -R g=u ./poetry.lock
chmod -R g=u ./data/cache

# Add current user to locusfocus group if not already part of group
if [[ -z $(groups $1 | grep -w locusfocus) ]]; then
    usermod -aG locusfocus $1
fi

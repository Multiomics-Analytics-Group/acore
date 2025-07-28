# #!/bin/bash
# Requires admin privileges to run for a GitHub organization
#
#
# Get all repos in the GitHub organization (replace ORGANIZATION with your org name)
repos=$(gh api /orgs/RasmussenLab/repos?sort=full_name --jq '.[].full_name')

# Calculate the ISO 8601 timestamp for 3 hours ago
ago=$(date -r $(($(date +%s) - 10800)) -Iseconds)

# Loop through each repository
# for repo in $repos; do
repo="RasmussenLab/pimms"
echo "Checking repo $repo"

# Get the check suite IDs for workflow runs created since $ago
check_suite_ids=$(gh api "/repos/$repo/actions/runs?created=>=$ago" --jq .workflow_runs.[].check_suite_id)

# Loop through each check suite ID
for check_suite_id in $check_suite_ids; do
# Get the annotation URLs for each check run in the suite
    annotation_urls=$(gh api "/repos/$repo/check-suites/$check_suite_id/check-runs" --jq .check_runs.[].output.annotations_url)

        # Loop through each annotation URL
        for url in $annotation_urls; do
            # Fetch and print each annotation message (warnings/errors)
            message=$(gh api "$url" --jq .[].message)
            if [[ -n "$message" ]]; then
            echo "- $message"
            fi
        done
done
echo ""
# done



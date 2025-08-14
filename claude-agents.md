# Cost-Efficient Claude Agents

## Quick Search Agent
Use when: Need to find specific code/files quickly
```
Search for X in *.py files, return only file paths
```

## Batch Fix Agent
Use when: Multiple similar fixes needed
```
Fix all occurrences of X with Y in [files]. No explanations.
```

## Code Review Agent
Use when: Need quick validation
```
Check [file] for: syntax errors, unused imports, type issues. List only problems.
```

## Git Commit Agent
Use when: Creating commits
```
Stage changes, commit with message describing what changed. No details.
```

## Debug Agent
Use when: Error occurs
```
Error: [paste error]. Find cause in [context]. Solution only.
```

## Refactor Agent
Use when: Optimizing code
```
Refactor [function/file] for: performance/readability. Changes only.
```

## Test Runner Agent
Use when: Running tests
```
Run tests, report failures only
```

## Container Check Agent
Use when: Checking container resources
```
List container sizes in cache. Show large ones (>100MB)
```

## Performance Agent
Use when: Checking workflow performance
```
Find slow processes in last run. Show top 3 with times.
```

## Dependency Agent
Use when: Managing dependencies
```
Check if [package] is used. Where? Can remove?
```

# Usage Tips for Cost Efficiency

1. **Be specific** - Exact file names, function names, error messages
2. **Request output only** - Skip explanations unless needed
3. **Batch operations** - Combine multiple similar tasks
4. **Use grep/glob first** - Before using general search
5. **Limit scope** - Specify exact directories/files
6. **Skip confirmations** - "Just do X" instead of "Can you do X?"
7. **One-line responses** - Request "answer only" or "result only"

# Examples

## Expensive (many tokens):
"Can you help me understand what this error means and walk me through fixing it?"

## Efficient (few tokens):
"Fix: [error]. File: main.py:45"

## Expensive:
"Please search through the codebase and find all Python files that might contain database connections"

## Efficient:
"Grep 'connect\|database\|db' in *.py"

## Expensive:
"Could you check if the container is cached and tell me about its size?"

## Efficient:
"ls -lh singularity_cache/*python*"
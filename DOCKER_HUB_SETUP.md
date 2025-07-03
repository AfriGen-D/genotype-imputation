# Docker Hub Setup for GitHub Actions

## Required Secrets

To enable pushing to Docker Hub, add these secrets to your GitHub repository:

### 1. Navigate to Repository Settings
- Go to your GitHub repository
- Click "Settings" tab
- Click "Secrets and variables" → "Actions"

### 2. Add Docker Hub Secrets

Add these repository secrets:

| Secret Name | Value | Description |
|-------------|-------|-------------|
| `DOCKERHUB_USERNAME` | Your Docker Hub username | Username for Docker Hub login |
| `DOCKERHUB_TOKEN` | Your Docker Hub access token | Access token (not password!) |

### 3. Create Docker Hub Access Token

1. Go to [Docker Hub](https://hub.docker.com/)
2. Login to your account
3. Click your username → "Account Settings"
4. Click "Security" tab
5. Click "New Access Token"
6. Name it "GitHub Actions" 
7. Select "Read, Write, Delete" permissions
8. Copy the token and add as `DOCKERHUB_TOKEN` secret

## Container Availability

Once set up, containers will be available on both registries:

### GitHub Container Registry (GHCR)
```bash
docker pull ghcr.io/afrigen-d/analysis:bcftools-1.20
docker pull ghcr.io/afrigen-d/plink2:plink2-2.0.0
```

### Docker Hub
```bash
docker pull afrigen-d/analysis:bcftools-1.20
docker pull afrigen-d/plink2:plink2-2.0.0
```

## Benefits of Dual Registry

- **Accessibility**: Docker Hub is more widely accessible
- **Redundancy**: Multiple sources for container images  
- **Discoverability**: Docker Hub has better search/discovery
- **Enterprise**: Some organizations prefer specific registries 
# https://github.com/vercel/next.js/blob/canary/examples/with-docker/Dockerfile
# e.g., docker build --build-arg NEXT_PUBLIC_GRAPHQL_ENDPOINT="http://localhost:9033" --build-arg NEXT_PUBLIC_PUBLIC_EMAIL_ADDRESS="foo@bar.com" -t hawk-eye-next next/

FROM node:22-alpine AS base

FROM base AS deps
# Check https://github.com/nodejs/docker-node/tree/b4117f9333da4138b03a546ec926ef50a31506c3#nodealpine to understand why libc6-compat might be needed.
RUN apk add --no-cache libc6-compat
WORKDIR /app

COPY ./locus-focus/package.json ./locus-focus/package-lock.json ./
RUN npm ci

# Rebuild the source code only when needed
FROM base AS builder
WORKDIR /app

ARG NEXT_PUBLIC_BROWSER_API_HOST
ARG NEXT_PUBLIC_APP_ENV

ENV NEXT_PUBLIC_BROWSER_API_HOST=${NEXT_PUBLIC_BROWSER_API_HOST}
ENV NEXT_PUBLIC_APP_ENV=${NEXT_PUBLIC_APP_ENV}

COPY --from=deps /app/node_modules ./node_modules
COPY ./locus-focus .

ENV NEXT_TELEMETRY_DISABLED=1

RUN npm run build

RUN ls -la

RUN cp -r .next/static .next/standalone/.next/

# Production image, copy all the files and run next
FROM base AS runner
WORKDIR /app

ARG NEXT_PORT
ENV PORT=${NEXT_PORT}

ENV NODE_ENV=production

ENV NEXT_TELEMETRY_DISABLED=1

RUN addgroup --system --gid 1001 nodejs
RUN adduser --system --uid 1001 nextjs

COPY --from=builder --chown=nextjs:nodejs /app/.next/standalone ./
COPY --from=builder --chown=nextjs:nodejs /app/.next/static ./.next/static

USER nextjs

ENV HOSTNAME="0.0.0.0"

CMD ["node", "server.js"]


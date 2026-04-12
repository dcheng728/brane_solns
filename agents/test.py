import asyncio
import os

# Remove invalid API key so claude-code-sdk uses claude.ai OAuth (Max subscription)
os.environ.pop("ANTHROPIC_API_KEY", None)

from claude_code_sdk import query, ClaudeCodeOptions, AssistantMessage, ResultMessage
from claude_code_sdk._errors import MessageParseError


async def main():
    try:
        async for message in query(
            prompt="What is the meaning of life",
            options=ClaudeCodeOptions(allowed_tools=["Read", "Edit", "Bash"]),
        ):
            if isinstance(message, AssistantMessage):
                for block in message.content:
                    if hasattr(block, "text"):
                        print(block.text)
            elif isinstance(message, ResultMessage):
                print(f"Done: {message.subtype}")
    except MessageParseError:
        pass  # Skip unknown message types (e.g. rate_limit_event)


asyncio.run(main())

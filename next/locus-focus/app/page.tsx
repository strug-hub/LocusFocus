import { Suspense } from "react";
import { ColocPage } from "@/components";

export default function Page() {
  const tissuesV8 = fetch(
    `${process.env.NEXT_PUBLIC_SERVER_API_HOST}/gtex/v8/tissues_list`
  ).then((r) => r.json());

  const tissuesV10 = fetch(
    `${process.env.NEXT_PUBLIC_SERVER_API_HOST}/gtex/v10/tissues_list`
  ).then((r) => r.json());

  return (
    <Suspense fallback={<div>Loading...</div>}>
      <ColocPage _tissuesV8={tissuesV8} _tissuesV10={tissuesV10} />
    </Suspense>
  );
}
